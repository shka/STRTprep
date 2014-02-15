require 'parallel'

tmp = Array.new
LIBIDS.each { |libid|
  taskid = "alignment_#{libid}"
  align = "out/stat/#{libid}.alignment.txt"
  spike = "out/stat/#{libid}.spike.txt.xz"
  task taskid => [align, spike]
  tmp.push(taskid)
  #
  tidxorg = File.expand_path(CONF[libid]['TRANSCRIPT']+'.1.ebwt')
  tidx = mytemppath('tidx')+'/'+tidxorg.pathmap('%f')
  file tidx => tidxorg do |t|
    newdir = t.name.pathmap('%d')
    sh "mkdir -p #{newdir}"
    sh "cp -r #{t.prerequisites[0].sub('.1.ebwt', '')}* #{newdir}"
  end
  gidxorg = File.expand_path(CONF[libid]['GENOMESPIKERIBO']+'.1.ebwt')
  gidx = mytemppath('gidx')+'/'+gidxorg.pathmap('%f')
  file gidx => gidxorg do |t|
    newdir = t.name.pathmap('%d')
    sh "mkdir -p #{newdir}"
    sh "cp -r #{t.prerequisites[0].sub('.1.ebwt', '')}* #{newdir}"
  end
  #
  # bams = Array.new
  seqs = Array.new
  open(File.expand_path(CONF[libid]['LAYOUT'])).each { |line|
    well, barcodegap = line.rstrip.split(/\t/)
    bam = "out/ali/#{libid}.#{well}/accepted_hits.bam"
    seq = "tmp/seq/#{libid}.#{well}.fq.gz"
    file bam => [seq, gidx, tidx] do |t|
      alignment(t)
    end
    seqs.push(seq)
    # bams.push(bam)
  }
  #
  # file align => bams do |t|
  #   stat_alignment(t)
  # end
  file align => seqs do |t|
    bams = Parallel.map(t.prerequisites, :in_threads => PROCS/2) { |seq|
      bam = seq.sub('tmp/seq', 'out/ali').sub('.fq.gz', '/accepted_hits.bam')
      Rake::Task[bam].invoke unless File.exist?(bam)
      bam
    }
    stat_alignment(t.name, bams)
  end
  #
  # file spike => bams do |t|
  #   stat_spike(t)
  # end
  file spike => seqs do |t|
    bams = Array.new
    seqs.each { |seq|
      bams.push(seq.sub('tmp/seq', 'out/ali').sub('.fq.gz', '/accepted_hits.bam'))
    }
    stat_spike(t.name, bams)
  end
}

task :alignment => tmp

def alignment(t)
  outdir = t.name.sub('/accepted_hits.bam', '')
  sh "mkdir -p #{outdir}"
  sh <<"PROCESS"
tophat --transcriptome-index #{t.prerequisites[2].sub('.1.ebwt', '')} \\
       --library-type fr-secondstrand --min-anchor 5 --coverage-search \\
       --output-dir #{outdir} --num-threads 2 --bowtie1 \\
       #{t.prerequisites[1].sub('.1.ebwt', '')} #{t.prerequisites[0]}
PROCESS
end

def stat_alignment(name, bams)
  libid = name.pathmap("%n").pathmap("%n")
  well2cnts = Hash.new
  Parallel.each(bams, :in_threads => PROCS) { |mbam|
    tmp = Array.new
    fp = open(mbam.sub('accepted_hits.bam', 'align_summary.txt'))
    fp.gets
    tmp.push(/Input\s*\:\s+(\d+)/.match(fp.gets).to_a[1])
    mapped = /Mapped\s*\:\s+(\d+)/.match(fp.gets).to_a[1]
    multi = /these\s*\:\s+(\d+)/.match(fp.gets).to_a[1]
    tmp.push(mapped.to_i-multi.to_i)
    tmp.push(mapped)
    fp.close
    open("| samtools view #{mbam} | cut -f 1,3,12- | grep '\tRNA_SPIKE_' | grep XS:A:+ | cut -f 1 | sort -u | wc -l").each { |line|
      tmp.push(line.rstrip) if line != "\n"
    }
    open("| samtools view #{mbam} | cut -f 1,3 | grep -E '\tRIBO_' | cut -f 1 | sort -u | wc -l").each { |line|
      tmp.push(line.rstrip) if line != "\n"
    }
    well = /#{libid}\.([^\/]+)/.match(mbam).to_a[1]
    well2cnts[well] = [libid, well] + tmp
  }
  # open("| samtools view #{mbam} | cut -f 1 | wc -l").each { |line|
  #   well2cnts[well].push(line.rstrip) if line != "\n"
  # }
  fp = open(name, 'w')
  fp.puts ['LIB', 'WELL', 'TOTAL', 'MAPPED.UNIQUE', 'MAPPED.ALL', 'SPIKE', 'RIBOSOMAL'].join("\t")
  well2cnts.keys.sort.each { |well| fp.puts well2cnts[well].join("\t") }
  fp.close
end

def stat_spike(name, bams)
  libid = name.pathmap("%n").pathmap("%n").pathmap("%n")
  fp = open("| xz -c --extreme > #{name}", 'w')
  bams.each { |bam|
    head = [libid, /#{libid}\.([^\/]+)/.match(bam).to_a[1], ''].join("\t")
    open("| samtools view #{bam} | cut -f 1,3,4,12- | grep '\tRNA_SPIKE_' | grep XS:A:+ | cut -f 2,3").each { |line| fp.puts head + line }
  }
  fp.close
end
