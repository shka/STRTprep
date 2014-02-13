require 'levenshtein'
require 'mkfifo'
require 'parallel'
require 'tempfile'
require 'yaml'

##

tmplibid = /STRTprep.([^\/.]+)/.match(Dir.pwd).to_a[1]
LIBID = ENV.key?('LIBID') ? ENV['LIBID'] : (tmplibid.nil? ? 'TMP' : tmplibid)

CONF = YAML.load_file('conf.yaml')
EBWT_PHYX = File.expand_path(CONF['INPUTS']['EBWTS']['PHYX'])
FASTQS = Array.new
CONF['INPUTS']['FASTQS'].each { |fastq| FASTQS.push(File.expand_path(fastq)) }

PROCS = ENV.key?('PROCS') ? ENV['PROCS'] : Parallel.processor_count

##

def read_barcodes
  wells = Array.new
  open('src/barcodes.txt').each do |line|
    wells.push(line.rstrip.split[0])
  end
  wells
end

WELLS = read_barcodes

##

task :default => [:alignment]

##

task :alignment => ['out/stat/alignment.txt', 'out/stat/spike.txt.gz']

BAMS = Array.new
WELLS.each { |well|
  bam = "out/ali/#{LIBID}.#{well}/accepted_hits.bam"
  file bam => ["out/seq/#{LIBID}.#{well}.fq.xz", 'src/knownGene.annotated.gtf', 'src/ebwt/ref.1.ebwt'] do |t|
    align(t)
  end
  BAMS.push(bam)
}

def align(t)
  outdir = t.name.sub('/accepted_hits.bam', '')
  logfile = t.name.sub('out/ali/', 'log/align.').sub('/accepted_hits.bam', '.log')
  gtf = Dir.exist?('tmp/knownGene.1.ebwt') ? '' : "-G #{t.prerequisites[1]}"
  sh "mkdir -p #{outdir}"
  sh <<"PROCESS"
xzcat #{t.prerequisites[0]} | \\
tophat #{gtf} --transcriptome-index tmp/knownGene --num-threads #{PROCS} \\
       --library-type fr-secondstrand --min-anchor 5 --coverage-search \\
       --bowtie1 --output-dir #{outdir} \\
       #{t.prerequisites[2].sub('.1.ebwt', '')} /dev/stdin > #{logfile} 2>&1
PROCESS
end

file 'out/stat/spike.txt.gz' => BAMS do |t|
  fp = open("| gzip -c --best > #{t.name}", 'w')
  t.prerequisites.each { |bam|
    head = [LIBID, /#{LIBID}\.([^\/]+)/.match(bam).to_a[1], ''].join("\t")
    open("| samtools view #{bam} | cut -f 1,3,4,12- | grep '\tRNA_SPIKE_' | grep XS:A:+ | cut -f 2,3").each { |line| fp.puts head + line }
  }
  fp.close
end

file 'out/stat/alignment.txt' => BAMS do |t|
  well2cnts = Hash.new
  Parallel.each(t.prerequisites, :in_threads => PROCS) { |mbam|
    tmp = Array.new
    fp = open(mbam.sub('accepted_hits.bam', 'align_summary.txt'))
    fp.gets
    tmp.push(/Input\:\s+(\d+)/.match(fp.gets).to_a[1])
    mapped = /Mapped\:\s+(\d+)/.match(fp.gets).to_a[1]
    multi = /these\:\s+(\d+)/.match(fp.gets).to_a[1]
    tmp.push(mapped.to_i-multi.to_i)
    tmp.push(mapped)
    fp.close
    open("| samtools view #{mbam} | cut -f 1,3,12- | grep '\tRNA_SPIKE_' | grep XS:A:+ | cut -f 1 | sort -u | wc -l").each { |line|
      tmp.push(line.rstrip) if line != "\n"
    }
    open("| samtools view #{mbam} | cut -f 1,3 | grep -E '\t(U13369|RIBOSOMAL)' | cut -f 1 | sort -u | wc -l").each { |line|
      tmp.push(line.rstrip) if line != "\n"
    }
    well = /#{LIBID}\.([^\/]+)/.match(mbam).to_a[1]
    well2cnts[well] = [LIBID, well] + tmp
  }
  # open("| samtools view #{mbam} | cut -f 1 | wc -l").each { |line|
  #   well2cnts[well].push(line.rstrip) if line != "\n"
  # }
  fp = open(t.name, 'w')
  fp.puts ['LIB', 'WELL', 'TOTAL', 'MAPPED.UNIQUE', 'MAPPED.ALL', 'SPIKE', 'RIBOSOMAL'].join("\t")
  well2cnts.keys.sort.each { |well| fp.puts well2cnts[well].join("\t") }
  fp.close
end

##

task :demultiplex => 'out/stat/demultiplex.txt'

WELLS.each { |well|
  file "out/seq/#{LIBID}.#{well}.fq.xz" => 'out/stat/demultiplex.txt'
}

file 'out/stat/demultiplex.txt' => ['src/barcodes.txt', "out/ali/#{LIBID}.phyX.bam"] do |t|

  sh 'mkdir -p out/seq out/stat'

  ofs =  4  # UMI length
  clen = 37 # cDNA length

  # barcode

  wells = Array.new
  bcs = Array.new
  bclens = Array.new
  open(t.prerequisites[0]).each do |line|
    cols = line.rstrip.split
    wells.push(cols[0])
    bcs.push(cols[1])
    bclens.push(cols[1].length)
  end

  bclens.uniq!
  if bclens.size != 1
    raise 'Barcode sequence lengths must be equivalent'
  end
  bclen = bclens[0]

  # non-phyX entries into multiple (fifo-)files

  fifo1paths = Array.new
  PROCS.times { |i| fifo1paths.push(mymkfifo('fifo1-')) }

  pid = Kernel.fork {
    fifo1s = Array.new
    fifo1paths.each { |fifo1path| fifo1s.push(open(fifo1path, 'w')) }
    total = 0
    open("| samtools view -f 4 -F 256 #{t.prerequisites[1]} | cut -f 1,10,11").each { |line|
      fifo1 = fifo1s[total % PROCS]
      fifo1.puts(line.rstrip)
      total += 1
      STDERR.write("Demultiplexing... #{total}\r") if total % 10000 == 0
    }
    fifo1s.each { |fifo1| fifo1.close }
    STDERR.puts
    Kernel.exit!
  }

  # barcode check

  fifo2paths = Array.new
  PROCS.times { |i|
    fifo2path = mymkfifo('fifo2-')
    fifo2paths.push(fifo2path)
    pid = Kernel.fork {
      open(fifo2path, 'w') { |fifo2|
        open(fifo1paths[i], 'r').each { |line|
          seqid, seq, qvs = line.rstrip.split(/\t/)
          tmpdists = Hash.new
          bcs.each_index { |bcidx|
            tmpdist = Levenshtein.distance(bcs[bcidx], seq[ofs, bclen], threshold=2)
            dist = tmpdist.nil? ? 2 : tmpdist
            tmpdists[bcidx] = dist
            break if dist < 2
          }
          dists = tmpdists.sort { |a, b| a[1] <=> b[1] }
          ## bc = dists[0][1] < 2 && dists[0][1] < dists[1][1] ? dists[0][0] : -1
          bc = dists[0][1] < 2 ? dists[0][0] : -1
          fifo2.puts([bc, seqid, seq, qvs].join("\t"))
        }
      }
      Kernel.exit!
    }
  }

  # [join] and write into files by barcode

  tmpwells = wells + ['nobc']

  fifo2s = Array.new
  fifo2paths.each { |fifo2path| fifo2s.push(open(fifo2path, 'r')) }
  fifo2done = Hash.new

  left = ofs+bclen
  right = clen > -1 ? ofs+bclen+clen-1 : -1
  outs = Array.new
  tmpwells.each_index { |i|
    well = tmpwells[i]
    if well == 'nobc'
      outs[i] = open("| xz -z -c -e > out/seq/#{LIBID}.nobc.fq.xz", 'w')
    else
      tmp = <<"DEDUPandFORMAT"
| sort -k 1 -r | cut -f 2- | uniq -f 2 \\
| ruby -F'\\t' -anle 'puts(["@"+$F[0], $F[1][#{left}..-1], "+", $F[2][#{left}..-1]].join("\\n"))' \\
| xz -z -c -e > out/seq/#{LIBID}.#{well}.fq.xz
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
    open("| xzcat out/seq/#{LIBID}.#{well}.fq.xz | grep ^@ | wc -l").each { |line|
      ucnts[well] = line.rstrip.to_i unless line.nil?
    }
  }

  open('out/stat/demultiplex.txt', 'w') { |fp|
    fp.puts ['LIB', 'WELL', 'BARCODE', 'TOTAL', 'UTOTAL'].join("\t")
    total = 0
    utotal = 0
    cnts.keys.sort.each { |well|
      fp.puts [LIBID, well, (well != 'nobc' ? bcs[wells.index(well)] : ''), cnts[well], ucnts[well]].join("\t")
      total += cnts[well]
      utotal += ucnts[well]
    }
    fp.puts ['TOTAL', '', '', total, utotal].join("\t")
  }

end

##

task :removePhyX => "out/ali/#{LIBID}.phyX.bam"

file "out/ali/#{LIBID}.phyX.bam" => [EBWT_PHYX+".1.ebwt"] + FASTQS do |t|
  sh 'mkdir -p out/ali log'
  tmp = Array.new
  sh "(gunzip -c #{t.prerequisites[1..-1].join(' ')} | bowtie -S -p #{PROCS} -t #{t.prerequisites[0].sub('.1.ebwt', '')} - | samtools view -@ #{PROCS} -S -b -o #{t.name} -) 2> log/removePhyX.log "
end

##

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
