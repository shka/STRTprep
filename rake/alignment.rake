require 'parallel'

tmp = Array.new
LIBIDS.each { |libid|
  taskid = "alignment_#{libid}"
  align = "out/stat/#{libid}.alignment.txt"
  annotation = "out/stat/#{libid}.annotation.txt"
  spike_raw = "out/stat/#{libid}.spike.raw.txt.xz"
  spike = "out/stat/#{libid}.spike.txt"
  demultiplex = "out/stat/#{libid}.demultiplex.txt"
  task taskid => [align, spike, annotation]
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
    bams.push("out/ali/#{libid}.#{well}/accepted_hits.bam")
  }
  #
  file align => [demultiplex, gidx, tidx] do |t|
    Parallel.map(bams, :in_threads => PROCS/2) { |bam|
      if !File.exist?(bam) || File.mtime(bam) < File.mtime(demultiplex)
        outdir = bam.sub('/accepted_hits.bam', '')
        sh "mkdir -p #{outdir}"
        sh <<"PROCESS"
tophat --transcriptome-index #{t.prerequisites[2].sub('.1.ebwt', '')} \\
       --library-type fr-secondstrand --min-anchor 5 --coverage-search \\
       --output-dir #{outdir} --num-threads 2 --bowtie1 \\
       #{t.prerequisites[1].sub('.1.ebwt', '')} \\
       #{bam.sub('out/ali', 'tmp/seq').sub('/accepted_hits.bam', '.fq.gz')}
PROCESS
      end
    }
    stat_alignment(t.name, bams)
  end
  #
  file spike => align do |t|
    stat_spike(t.name, bams)
  end
  #
  file annotation => align do |t|
    tmp = Parallel.map(bams, :in_threads => PROCS) { |bam|
      beds = Array.new
      0.upto(8) { |i|
        beds.push(bam.sub('out/', 'tmp/').sub('/accepted_hits.bam', "/nonclass#{i}.bed.gz"))
      }
      Rake::Task[beds[-1]].invoke
      cnts = Array.new
      0.upto(8) { |i| cnts.push(`gunzip -c #{beds[i]} | wc -l`.strip.to_i) }
      ret = [bam.pathmap('%3d').pathmap('%f').sub('.', "\t"), cnts[0]]
      1.upto(8) { |i| ret.push(cnts[i-1]-cnts[i]) }
      ret.push(ret[-9]-ret[-1]-ret[-2]-ret[-3]-ret[-4]-ret[-5]-ret[-6]-ret[-7]-ret[-8])
      ret
    }
    fp = open(t.name, 'w')
    fp.puts ['LIB', 'WELL', 'MAPPED.UNIQ.NONSPIKE.NONRIBO', "5'-UTR", 'CODING.UPSTREAM', 'CDS', "3'-UTR", 'NONCODING.1ST', 'NONCODING.UPSTREAM', 'NONCODING.OTHER', 'INTRON', 'UNANNOTATED'].join("\t")
    tmp.each { |row| fp.puts row.join("\t") }
    fp.close
  end
}

task :alignment => tmp

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

rule /tmp\/ali\/.*\/nonclass0\.bed\.gz/ => proc { |t|
  t.sub('tmp/', 'out/').sub('/nonclass0.bed.gz', '/accepted_hits.bam')
} do |t|
  sh "mkdir -p #{t.name.pathmap('%d')}"
  bam = t.prerequisites[0]
  #
  acc2dups = Hash.new
  open("| samtools view #{bam} | cut -f 1 | sort | uniq -c").each { |line|
    dups, acc = line.strip.split(/\s/)
    acc2dups[acc] = '' if dups != '1'
  }
  #
  out = open("| gzip -c --best > #{t.name}", 'w')
  open("| bamToBed -i #{bam}").each { |line|
    cols = line.rstrip.split(/\t/)
    next if acc2dups.key?(cols[3]) || (cols[0] =~ /^RNA_SPIKE_/ && cols[5] == '+') || cols[0] =~ /^RIBO_/
    if cols[5] == '+'
      out.puts([cols[0], cols[1], cols[1].to_i+1, cols[3], 1, '+'].join("\t"))
    else
      out.puts([cols[0], cols[2].to_i-1, cols[2], cols[3], 1, '-'].join("\t"))
    end
  }
  out.close
end

rule /tmp\/ali\/.*\/nonclass[1-8]\.bed\.gz/ => proc { |t|
  libid = t.pathmap('%3d').pathmap('%n')
  clsid = /nonclass(\d)/.match(t).to_a[1].to_i
  [t.sub("nonclass#{clsid}", "nonclass#{clsid-1}"),
    File.expand_path(t.sub(t.pathmap('%3d')+'/non',
                           CONF[libid]['TRANSCRIPT']+'.'))]
} do |t|
  sh "bedtools intersect -s -a #{t.prerequisites[0]} -b #{t.prerequisites[1]} -v | gzip -c > #{t.name}"
end
