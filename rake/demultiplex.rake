require 'levenshtein'
require 'mkfifo'
require 'parallel'
require 'tempfile'

tmp = Array.new
LIBIDS.each { |libid|
  taskid = "demultiplex_#{libid}"
  tmp.push(taskid)
  task taskid => "out/stat/#{libid}.demultiplex.txt"
  # WELLS.each { |well|
  #   file "out/seq/#{libid}.#{well}.fq.xz" => 'out/stat/demultiplex.txt'
  # }
}

task 'demultiplex' => tmp

rule '.demultiplex.txt' => proc { |target|
  libid = target.pathmap("%n").pathmap("%n")
  [File.expand_path(CONF[libid]['LAYOUT']),
   target.sub("out/stat/#{libid}.demultiplex.txt", "out/ali/#{libid}.phyX.bam")]
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
      open(fifo2path, 'w') { |fifo2|
        open(fifo1paths[i], 'r').each { |line|
          seqid, seq, qvs = line.rstrip.split(/\t/)
          tmpdists = Hash.new
          barcodegaps.each_index { |bcidx|
            tmpdist = Levenshtein.distance(barcodegaps[bcidx], seq[umi, barcode+gap], threshold=2)
            dist = tmpdist.nil? ? 2 : tmpdist
            tmpdists[bcidx] = dist
            break if dist < 2
          }
          dists = tmpdists.sort { |a, b| a[1] <=> b[1] }
          bc = dists[0][1] < 2 ? dists[0][0] : -1
          fifo2.puts([bc, seqid, seq, qvs].join("\t"))
        }
      }
      Kernel.exit!
    }
  }

  # [join] and write into files by barcode

  tmpwells = wells + ['nonbc']

  fifo2s = Array.new
  fifo2paths.each { |fifo2path| fifo2s.push(open(fifo2path, 'r')) }
  fifo2done = Hash.new

  left = umi+barcode+gap
  right = cdna > -1 ? umi+barcode+gap+cdna-1 : -1
  outs = Array.new
  tmpwells.each_index { |i|
    well = tmpwells[i]
    if well == 'nonbc'
      outs[i] = open("| xz -z -c -e > out/seq/#{libid}.#{well}.fq.xz", 'w')
    else
      tmp = <<"DEDUPandFORMAT"
| sort -k 1 -r | cut -f 2- | uniq -f 2 \\
| ruby -F'\\t' -anle 'puts(["@"+$F[0], $F[1][#{left}..-1], "+", $F[2][#{left}..-1]].join("\\n"))' \\
| gzip -c > tmp/seq/#{libid}.#{well}.fq.gz
DEDUPandFORMAT
      outs[i] = open(tmp, 'w')
    end
  }

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
        if bc == -1
          outs[-1].puts(['@'+seqid, seq, '+', qvs].join("\n"))
        else
          f1 = seq[0..right]
          f2 = qvs[0..right]
          outs[bc].puts([f1+f2, seqid, f1, f2].join("\t"))
        end
      end
    end
    if fifo2done.size == fifo2s.size
      break
    end
  }
  fifo2s.each { |fifo2| fifo2.close }
  outs.each { |out| out.close }

  Process.waitall

  ucnts = Hash.new
  Parallel.each(cnts.keys, :in_threads => PROCS) { |well|
    ucnts[well] = 0
    if well == 'nonbc'
      ucnts[well] = cnts[well]
    else
      open("| zcat tmp/seq/#{libid}.#{well}.fq.gz | wc -l").each { |line|
        ucnts[well] = line.rstrip.to_i/4 unless line.nil?
      }
    end
  }

  open("out/stat/#{libid}.demultiplex.txt", 'w') { |fp|
    fp.puts ['LIB', 'WELL', 'BARCODE', 'TOTAL', 'UTOTAL'].join("\t")
    total = 0
    utotal = 0
    cnts.keys.sort.each { |well|
      fp.puts [libid, well, (well != 'nonbc' ? barcodegaps[wells.index(well)] : ''), cnts[well], ucnts[well]].join("\t")
      total += cnts[well]
      utotal += ucnts[well]
    }
    fp.puts ['TOTAL', '', '', total, utotal].join("\t")
  }
end

#

@mytemppaths = Array.new

def mytemppath(basename, tmpdir = Dir::tmpdir)
  fp = Tempfile.open(basename, tmpdir)
  path = fp.path
  @mytemppaths.push(path)
  fp.close!
  path
end

END { @mytemppaths.each { |path| File.unlink(path) if File.exist?(path) } }

##

def mymkfifo(basename, tmpdir = Dir::tmpdir)
  path = mytemppath(basename, tmpdir)
  File.mkfifo(path)
  path
end
