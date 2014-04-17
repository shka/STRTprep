require 'parallel'

####
#
# dummy paths
#
tmp = Array.new
gidxmap = Hash.new
tidxmap = Hash.new
LIBIDS.each { |libid|
  bams = Array.new
  beds4annotation = Array.new
  open(File.expand_path(CONF[libid]['LAYOUT'])).each { |line|
    cols = line.rstrip.split
    wellid = cols[0]
    bams.push("out/ali/#{libid}.#{wellid}/accepted_hits.bam")
    beds4annotation.push("tmp/ali/#{libid}.#{wellid}/nonclass8.bed.gz")
  }
  #
  timestamp = "tmp/#{libid}.alignment.timestamp"
  file timestamp => beds4annotation do |t|
    sh "touch #{t.name}"
    begin
      Parallel.each(t.prerequisites, :in_threads => PROCS) { |target|
        Rake::Task[target].invoke
      }
    rescue
      sh "rm -rf #{timestamp}"
    end
  end
  #
  gidxorg = File.expand_path(CONF[libid]['GENOMESPIKERIBO']+'.1.ebwt')
  gidxdir = mytemppath('gidx')
  gidx = gidxdir+'/'+gidxorg.pathmap('%f')
  sh "mkdir -p #{gidxdir}; cp -rp #{gidxorg.sub('.1.ebwt', '')}* #{gidxdir}"
  gidxmap[gidxorg] = gidx
  #
  tidxorg = File.expand_path(CONF[libid]['TRANSCRIPT']+'.1.ebwt')
  tidxdir = mytemppath('tidx')
  tidx = tidxdir+'/'+tidxorg.pathmap('%f')
  sh "mkdir -p #{tidxdir}; cp -rp #{tidxorg.sub('.1.ebwt', '')}* #{tidxdir}"
  tidxmap[tidxorg] = tidx
  #
  stat_alignment = "out/stat/#{libid}.alignment.txt"
  file stat_alignment => [timestamp] + bams do |t|
    report_alignment(t.name, t.prerequisites[1..-1])
  end
  tmp.push(stat_alignment)
  #
  stat_spike = "out/stat/#{libid}.spike.txt"
  file stat_spike => [timestamp] + bams do |t|
    report_spike(t.name, t.prerequisites[1..-1])
  end
  tmp.push(stat_spike)
  #
  stat_annotation = "out/stat/#{libid}.annotation.txt"
  file stat_annotation => [timestamp] + bams do |t|
    report_annotation(t.name, t.prerequisites[1..-1])
  end
  tmp.push(stat_annotation)
}
task :alignment => tmp + [:demultiplex]

####
#
# file out/ali/#{libid}.#{wellid}/accepted_hits.bam
#      => [tmp/seq/#{libid}.#{wellid}.fq.gz,
#          #{libid.GENOMESPIKERIBO}.1.ebwt,
#          #{libid.TRANSCRIPT}.1.ebwt]
#
rule /out\/ali\/[^\/]+\/accepted_hits\.bam/ => proc { |target|
  libid, wellid = /\/([^.\/]+)\.([^.\/]+)\//.match(target)[1..2]
  [ "tmp/seq/#{libid}.#{wellid}.fq.gz",
    File.expand_path(CONF[libid]['GENOMESPIKERIBO']+'.1.ebwt'),
    File.expand_path(CONF[libid]['TRANSCRIPT']+'.1.ebwt')]
} do |t|
  outdir = t.name.sub('/accepted_hits.bam', '')
  sh "mkdir -p #{outdir}"
  sh <<"PROCESS"
tophat \\
  --transcriptome-index #{tidxmap[t.prerequisites[2]].sub('.1.ebwt', '')} \\
  --library-type fr-secondstrand --min-anchor 5 --coverage-search \\
  --output-dir #{outdir} --num-threads 2 --bowtie1 \\
    #{gidxmap[t.prerequisites[1]].sub('.1.ebwt', '')} #{t.prerequisites[0]}
PROCESS
end

####
#
def report_annotation(name, bams)
  tmp = Parallel.map(bams, :in_threads => PROCS) { |bam|
    beds = Array.new
    0.upto(8) { |i|
      beds.push(bam.sub('out/', 'tmp/').sub('/accepted_hits.bam', "/nonclass#{i}.bed.gz"))
    }
    cnts = Array.new
    0.upto(8) { |i| cnts.push(`gunzip -c #{beds[i]} | wc -l`.strip.to_i) }
    ret = [bam.pathmap('%3d').pathmap('%f').sub('.', "\t"), cnts[0]]
    1.upto(8) { |i| ret.push(cnts[i-1]-cnts[i]) }
    ret.push(ret[-9]-ret[-1]-ret[-2]-ret[-3]-ret[-4]-ret[-5]-ret[-6]-ret[-7]-ret[-8])
    ret
  }
  fp = open(name, 'w')
  fp.puts ['LIB', 'WELL', 'MAPPED.UNIQ.NONSPIKE.NONRIBO', "5'-UTR", 'CODING.UPSTREAM', 'CDS', "3'-UTR", 'NONCODING.1ST', 'NONCODING.UPSTREAM', 'NONCODING.OTHER', 'INTRON', 'UNANNOTATED', "5'-END.RATE"].join("\t")
  tmp.each { |row| fp.puts (row + [(row[-9]+row[-8]+row[-5]+row[-4]).to_f/row[-10].to_f]).join("\t") }
  fp.close
end

####
#
def report_alignment(name, bams)
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
    tmp.push(tmp[-1].to_f/tmp[-2].to_f)
    tmp.push(mapped)
    fp.close
    open("| samtools view #{mbam} | cut -f 1,3,12- | grep '\tRNA_SPIKE_' | grep XS:A:+ | cut -f 1 | sort -u | wc -l").each { |line|
      tmp.push(line.rstrip) if line != "\n"
    }
    open("| samtools view #{mbam} | cut -f 1,3 | grep -E '\tRIBO_' | cut -f 1 | sort -u | wc -l").each { |line|
      tmp.push(line.rstrip) if line != "\n"
    }
    tmp.push((tmp[-3].to_f-tmp[-2].to_f-tmp[-1].to_f)/tmp[-2].to_i)
    well = /#{libid}\.([^\/]+)/.match(mbam).to_a[1]
    well2cnts[well] = [libid, well] + tmp
  }
  fp = open(name, 'w')
  fp.puts ['LIB', 'WELL', 'TOTAL', 'MAPPED.UNIQUE', 'MAPPED.UNIQUE/TOTAL', 'MAPPED.ALL', 'SPIKE', 'RIBOSOMAL', '(MAPPED.ALL-SPIKE-RIBOSOMAL)/SPIKE'].join("\t")
  well2cnts.keys.sort.each { |well| fp.puts well2cnts[well].join("\t") }
  fp.close
end

####
#
def report_spike(name, bams)
  libid = name.pathmap("%n").pathmap("%n").pathmap("%n")
  raw = name.sub(/\.txt$/, '.raw.txt.xz')
  fp = open("| xz -c --extreme > #{raw}", 'w')
  Parallel.map(bams, :in_threads => PROCS) { |bam|
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
  lwsr2cnt.keys.sort.each { |libwell|
    sr2cnt = lwsr2cnt[libwell]
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

####
#
rule /tmp\/ali\/.*\/nonclass0\.bed\.gz/ => proc { |t|
  t.sub('tmp/', 'out/').sub('/nonclass0.bed.gz', '/accepted_hits.bam')
} do |t|
  sh "mkdir -p #{t.name.pathmap('%d')}"
  sh "bin/_process_nonclass0.rb #{t.prerequisites[0]} | gzip -c > #{t.name}"
end

####
#
rule /tmp\/ali\/.*\/nonclass[1-8]\.bed\.gz/ => proc { |t|
  libid = t.pathmap('%3d').pathmap('%n')
  clsid = /nonclass(\d)/.match(t).to_a[1].to_i
  [t.sub("nonclass#{clsid}", "nonclass#{clsid-1}"),
    File.expand_path(t.sub(t.pathmap('%3d')+'/non',
                           CONF[libid]['TRANSCRIPT']+'.'))]
} do |t|
  sh "bedtools intersect -s -a #{t.prerequisites[0]} -b #{t.prerequisites[1]} -v | gzip -c > #{t.name}"
end

####
#
# cleaning
#
task 'clean_alignment' do
  LIBIDS.each { |libid|
    sh "rm -rf tmp/ali/#{libid}.* out/ali/#{libid}.*/* out/stat/#{libid}.alignment.txt out/stat/#{libid}.annotation.txt out/stat/#{libid}.spike.* tmp/#{libid}.alignment.timestamp"
  }
end
