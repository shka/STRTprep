require 'levenshtein'
require 'parallel'

####
#
# dummy paths
#
tmp = Array.new
LIBIDS.each { |libid|
  timestamp = "tmp/#{libid}.demultiplex.timestamp"
  file "out/stat/#{libid}.demultiplex.txt"
  open(File.expand_path(CONF[libid]['LAYOUT'])).each { |line|
    cols = line.rstrip.split
    file "tmp/seq/#{libid}.#{cols[0]}.fq.gz" => timestamp
  }
  tmp.push(timestamp)
}
task :demultiplex => tmp

####
#
# file (tmp/#{libid}).demultiplex.timestamp
#      => [ #{libid.LAYOUT}, out/ali/#{libid}.phyX.bam ]
#
rule '.demultiplex.timestamp' => proc { |target|
  libid = target.pathmap("%n").pathmap("%n")
  [File.expand_path(CONF[libid]['LAYOUT']), "out/ali/#{libid}.phyX.bam"]
} do |t|

  sh 'mkdir -p out/seq out/stat tmp/seq'

  libid = t.name.pathmap("%n").pathmap("%n")
  conf = CONF[libid]

  umi = conf['UMI'].to_i
  barcode = conf['BARCODE'].to_i
  gap = conf['GAP'].to_i
  cdna = conf['CDNA'].to_i

  # layout

  wells = Array.new
  barcodegaps = Array.new
  open(t.prerequisites[0]).each do |line|
    cols = line.rstrip.split
    raise "Invalid barcode/gap length" if cols[1].length != barcode+gap
    wells.push(cols[0])
    barcodegaps.push(cols[1])
  end

  sh "touch tmp/#{libid}.demultiplex.timestamp"
  begin

    # [map] non-phyX entries into multiple (fifo-)files

    fifo1paths = Array.new
    PROCS.times { |i| fifo1paths.push(mymkfifo('fifo1-')) }

    pid = Kernel.fork {
      fifo1s = Array.new
      fifo1paths.each { |fifo1path| fifo1s.push(open(fifo1path, 'w')) }
      i = 0
      open("| samtools view -f 4 -F 256 #{t.prerequisites[1]} | cut -f 1,10,11").each { |line|
        fifo1s[i % PROCS].puts(line.rstrip)
        i += 1
      }
      fifo1s.each { |fifo1| fifo1.close }
      Kernel.exit!
    }

    # [mapped] barcode check

    fifo2paths = Array.new
    PROCS.times { |i|
      fifo2path = mymkfifo('fifo2-')
      fifo2paths.push(fifo2path)
      pid = Kernel.fork {
        qvcheck = Regexp.new("[!-#]") # QV2 == "#",  QV17 == "2"
        open(fifo2path, 'w') { |fifo2|
          open(fifo1paths[i], 'r').each { |line|
            seqid, seq, qvs = line.rstrip.split(/\t/)
            if qvcheck.match(qvs).nil?
              tmpdists = Hash.new
              barcodegaps.each_index { |bcidx|
                tmpdist = Levenshtein.distance(barcodegaps[bcidx],
                                               seq[umi, barcode+gap],
                                               threshold=2)
                dist = tmpdist.nil? ? 2 : tmpdist
                tmpdists[bcidx] = dist
                break if dist < 2
              }
              dists = tmpdists.sort { |a, b| a[1] <=> b[1] }
              bc = dists[0][1] < 2 ? dists[0][0] : -1
              fifo2.puts([bc, seqid, seq, qvs].join("\t"))
            else
              fifo2.puts([-2, seqid, seq, qvs].join("\t"))
            end
          }
        }
        Kernel.exit!
      }
    }

    # [join] and write into files by barcode

    tmpwells = wells + ['lowqv', 'nonbc']

    fifo2s = Array.new
    fifo2paths.each { |fifo2path| fifo2s.push(open(fifo2path, 'r')) }
    fifo2done = Hash.new

    left = umi+barcode+gap
    right = cdna > -1 ? umi+barcode+gap+cdna-1 : -1
    outs = Array.new
    tmpwells.each_index { |i|
      well = tmpwells[i]
      if well == 'nonbc' || well == 'lowqv'
        outs[i] = open("| gzip --best -c > out/seq/#{libid}.#{well}.fq.gz", 'w')
      else
        tmp = <<"DEDUPandFORMAT"
| sort -d -k 2,3 -r -S #{sprintf("%d%%", 100/PROCS)} \\
| ruby -F'\\t' -anle 'BEGIN { pre = ""}; if pre != $F[1] then puts(["@"+$F[0], $F[1][#{left}..-1], "+", $F[2][#{left}..-1]].join("\\n")); pre = $F[1] end' \\
| gzip -c > tmp/seq/#{libid}.#{well}.fq.gz
DEDUPandFORMAT
        outs[i] = open(tmp, 'w')
      end
    }

    # counting

    cnts = Hash.new
    fifo2s.cycle { |fifo2|
      unless fifo2done.key?(fifo2)
        line = fifo2.gets
        if line.nil?
          fifo2done[fifo2] = ''
        else
          bcstr, seqid, seq, qvs = line.rstrip.split(/\t/)
          bc = bcstr.to_i
          well = tmpwells[bc]
          cnts[well] = 0 unless cnts.key?(well)
          cnts[well] += 1
          if bc == -1 || bc == -2
            outs[bc].puts(['@'+seqid, seq, '+', qvs].join("\n"))
          else
            f1 = seq[0..right]
            f2 = qvs[0..right]
            outs[bc].puts([seqid, f1, f2].join("\t"))
          end
        end
      end
      if fifo2done.size == fifo2s.size
        break
      end
    }
    Parallel.each(fifo2s, :in_threads => PROCS) { |fifo2| fifo2.close }
    Parallel.each(outs, :in_threads => PROCS) { |out| out.close }

    Process.waitall

    # report

    ucnts = Hash.new
    Parallel.each(cnts.keys, :in_threads => PROCS) { |well|
      ucnts[well] = 0
      if well == 'nonbc' || well == 'lowqv'
        ucnts[well] = cnts[well]
      else
        open("| gunzip -c tmp/seq/#{libid}.#{well}.fq.gz | wc -l").each { |line|
          ucnts[well] = line.rstrip.to_i/4 unless line.nil?
        }
      end
    }

    open("out/stat/#{libid}.demultiplex.txt", 'w') { |fp|
      fp.puts ['LIB', 'WELL', 'BARCODE', 'TOTAL', 'UTOTAL', 'REDUNDANCY'].join("\t")
      total = 0
      utotal = 0
      cnts.keys.sort.each { |well|
        if well != 'nonbc' && well != 'lowqv'
          fp.puts [libid, well, barcodegaps[wells.index(well)], cnts[well], ucnts[well], cnts[well].to_f/ucnts[well].to_f].join("\t")
        else
          fp.puts [libid, well, '', cnts[well], ucnts[well]].join("\t")
        end
        total += cnts[well]
        utotal += ucnts[well]
      }
      fp.puts [libid, 'TOTAL', '', total, utotal].join("\t")
    }

  rescue
    sh "rm -rf tmp/#{libid}.demultiplex.timestamp"
  end
end

####
#
# cleaning
#
task 'clean_demultiplex' do
  LIBIDS.each { |libid|
    sh "rm -rf tmp/seq/#{libid}.*.fq.gz out/seq/#{libid}.*.fq.gz out/stat/#{libid}.demultiplex.txt tmp/#{libid}.demultiplex.timestamp"
  }
end
