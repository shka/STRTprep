require 'parallel'

tmp = Array.new
LIBIDS.each { |libid|
  taskid = "alignment_#{libid}"
  align = "out/stat/#{libid}.alignment.txt"
  spike_raw = "out/stat/#{libid}.spike.raw.txt.xz"
  spike = "out/stat/#{libid}.spike.txt"
  demultiplex = "out/stat/#{libid}.demultiplex.txt"
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
  bams = Array.new
  open(File.expand_path(CONF[libid]['LAYOUT'])).each { |line|
    well, barcodegap = line.rstrip.split(/\t/)
    bam = "out/ali/#{libid}.#{well}/accepted_hits.bam"
    file bam => [demultiplex, gidx, tidx] do |t|
      alignment(t)
    end
    bams.push(bam)
  }
  #
  file align do |t|
    Parallel.map(bams, :in_threads => PROCS/2) { |bam|
      Rake::Task[bam].invoke unless File.exist?(bam)
    }
    stat_alignment(t.name, bams)
  end
  #
  file spike => align do |t|
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
       #{t.prerequisites[1].sub('.1.ebwt', '')} \\
       #{t.name.sub('out/ali', 'tmp/seq').sub('/accepted_hits.bam', '.fq.gz')}
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
  fp = open(name, 'w')
  fp.puts ['LIB', 'WELL', 'TOTAL', 'MAPPED.UNIQUE', 'MAPPED.ALL', 'SPIKE', 'RIBOSOMAL'].join("\t")
  well2cnts.keys.sort.each { |well| fp.puts well2cnts[well].join("\t") }
  fp.close
end

def stat_spike(name, bams)
  libid = name.pathmap("%n").pathmap("%n").pathmap("%n")
  raw = name.sub(/\.txt$/, '.raw.txt.xz')
  fp = open("| xz -c --extreme > #{raw}", 'w')
  Parallel.map(bams, :in_threads => PROCS/2) { |bam|
    head = [libid, /#{libid}\.([^\/]+)/.match(bam).to_a[1], ''].join("\t")
    open("| samtools view #{bam} | cut -f 1,3,4,12- | grep '\tRNA_SPIKE_' | grep XS:A:+ | cut -f 2,3").each { |line| fp.puts head + line }
  }
  fp.close
  #
  lwsr2cnt = Hash.new
  s2cnt = Hash.new
  open("| xzcat #{raw}").each { |line|
    lib, well, spike, pos = line.rstrip.split(/\t/)
    libwell = "#{lib}\t#{well}"
    lwsr2cnt[libwell] = Hash.new unless lwsr2cnt.key?(libwell)
    lwsr2cnt[libwell][spike] = [0, 0, 0] unless lwsr2cnt[libwell].key?(spike)
    range = pos.to_i<=10 ? 0 : (pos.to_i<=100 ? 1 : 2)
    lwsr2cnt[libwell][spike][range] += 1
    if s2cnt.key?(spike)
      s2cnt[spike] += 1
    else
      s2cnt[spike] = 1
    end
  }
  spike_ordered = s2cnt.keys.sort
  fp = open(name, 'w')
  tmp = ['LIB', 'WELL', 'ORDER']
  spike_ordered.each { |spike|
    tmp.push("#{spike}, pos<=10", "#{spike}, pos<=100", "#{spike}, 100<pos", "#{spike}, pos<=10 rate")
  }
  fp.puts tmp.join("\t")
  lwsr2cnt.each { |libwell, sr2cnt|
    tmp = [libwell]
    tmp2 = sr2cnt.keys.sort { |a, b|
      r2cnt_a = sr2cnt[a]
      r2cnt_b = sr2cnt[b]
      r2cnt_b[0]+r2cnt_b[1]+r2cnt_b[2] <=> r2cnt_a[0]+r2cnt_a[1]+r2cnt_a[2]
    }
    tmp.push(tmp2.join(','))
    spike_ordered.each { |spike|
      if sr2cnt.key?(spike)
        r2cnt = sr2cnt[spike]
        tmp = tmp + r2cnt + [r2cnt[0].to_f/(r2cnt[0]+r2cnt[1]+r2cnt[2])]
      else
        tmp.push(0, 0, 0, 0)
      end
    }
    fp.puts tmp.join("\t")
  }
  fp.close
end
