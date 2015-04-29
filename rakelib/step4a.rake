##
## Step 4a - transcript assembly & the annotation
##

require 'csv'

def step4a_transcripts_sources(t)
  cls = t.pathmap("%d").pathmap("%f").pathmap("%X")
  src = ['src/samples.csv']
  samples = CSV.table(src[0])
  samples.each do |row|
    if row[:classtfe] == cls
      src.push("out/bam/#{row[:library]}.#{row[:well]}.bam")
    end
  end
  return src
end

rule /\.step4a\/transcripts\.gtf$/ => [->(path){ step4a_transcripts_sources(path) }] do |t|
  dir = t.name.pathmap("%d")
  mkdir_p dir
  bam = "#{dir}/merged.bam"
  bams = t.sources[1..-1].sort.join(' ')
  sh "samtools merge -f #{bam} #{bams} > #{t.name}.log 2>&1"
  tmp = "#{dir}/merged.bam.bak"
  if (!File.exist?(t.name) ||
      !File.exist?(tmp) ||
      `gmd5sum #{bam} #{tmp} | gcut -d ' ' -f 1 | guniq | gwc -l`.to_i != 1)
    cls = t.name.pathmap("%d").pathmap("%f").pathmap("%X")
    sh "(cufflinks -o #{dir} -p #{PROCS} --library-type fr-secondstrand -L #{cls} #{bam}) >> #{t.name}.log 2>&1"
    sh "cp -p #{bam} #{tmp}"
  else
    puts "... skipped #{t.name}, since the qualified samples were identical with the previous run"
    sh "touch #{t.name}"
  end
end

##

rule /\.step4a\/firstExons\.bed\.gz$/ => [
  ->(path){ return path.sub(/\/firstExons\.bed\.gz$/, '/transcripts.gtf') }] do |t|
  outfp = open("| gzip -c > #{t.name}", 'w')
  infp = open(t.source)
  tid = nil
  str = nil
  buf = nil
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    if cols[2] == 'transcript'
      outfp.puts buf unless buf.nil?
      tid = /transcript_id "([^"]+)"/.match(cols[8]).to_a[1]
      str = cols[6]
      buf = nil
    elsif cols[2] == 'exon'
      if buf.nil? || str == '-'
        buf = [cols[0], cols[3].to_i-1, cols[4], tid, 0, str].join("\t")
      end
    end
  end
  infp.close
  outfp.puts buf unless buf.nil?
  outfp.close
end

##

def step4a_fivePrimes_sources(t)
  cls = t.pathmap("%d").pathmap("%f").pathmap("%X")
  src = Array.new
  samples = CSV.table('src/samples.csv')
  samples.each do |row|
    if row[:classtfe] == cls
      src.push("tmp/#{row[:library]}.#{row[:well]}.step2g")
    end
  end
  return src
end

rule /\.step4a\/fivePrimes.bed.gz$/ => [->(path){ step4a_fivePrimes_sources(path) }] do |t|
  sh <<EOF
gunzip -c #{t.sources.join(' ')} \
| gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -t '\t' -k 1,1 -k 2,2n \
| mergeBed -s -c 4,5,6 -o distinct,mean,distinct -d -1 -i - \
| gzip -c > #{t.name}
EOF
end
  
##

step4a_firstExons = Array.new
step4a_fivePrimes = Array.new
begin 
  samples = CSV.table('src/samples.csv')
  classes = Array.new
  samples.each { |row| classes.push(row[:classtfe]) if row[:classtfe] != 'NA' }
  classes.uniq.each do |cls|
    step4a_firstExons.push("tmp/#{cls}.step4a/firstExons.bed.gz")
    step4a_fivePrimes.push("tmp/#{cls}.step4a/fivePrimes.bed.gz")
  end
end

##

file 'out/byTFE/regions.bed.gz' => step4a_firstExons do |t|
  mkdir_p t.name.pathmap('%d')
  
  outfp = open("| gzip -c > #{t.name}", 'w')
  infp = open("| gunzip -c #{t.sources.join(' ')} | gsort -t '\t' -k 1,1 -k 2,2n | mergeBed -s -c 6 -o distinct -i -")
  spikes = Hash.new
  tfe = 0
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    if /^RNA_SPIKE_/ =~ cols[0] && cols[3] == '+' && !spikes.key?(cols[0])
      outfp.puts [cols[0], cols[1], cols[2], cols[0], 0, cols[3]].join("\t")
      spikes[cols[0]] = ''
    else
      outfp.puts [cols[0], cols[1], cols[2], "TFE#{tfe}", 0, cols[3]].join("\t")
      tfe = tfe+1
    end
  end
  infp.close
  outfp.close
end

##

file 'tmp/byTFE/fivePrimes.bed.gz' => step4a_fivePrimes do |t|
  sh <<EOF
gunzip -c #{t.sources.join(' ')} \
| gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -t '\t' -k 1,1 -k 2,2n \
| mergeBed -s -c 4,5,6 -o distinct,mean,distinct -d -1 -i - \
| gzip -c > #{t.name}
EOF
end

##

file 'out/byTFE/peaks.bed.gz' => ['out/byTFE/regions.bed.gz',
                                  'tmp/byTFE/fivePrimes.bed.gz'] do |t|
  sh <<EOF
intersectBed -wa -wb -s -a #{t.source} -b #{t.sources[1]} \
| gcut -f 4,7,8,9,11,12 \
| gawk 'BEGIN{ FS="\t"; OFS="\t" }{p=$6=="+"?$3:-$4;print $2,$3,$4,$1,$5,$6,p,$1}' \
| gsort --parallel=#{PROCS} -S #{50/(PROCS+1)}% -t '\t' -k 8,8 -k 5,5gr -k 7,7g \
| guniq -f 7 \
| gcut -f 1-6 \
| gsort --parallel=#{PROCS} -S #{50/(PROCS+1)}% -t '\t' -k 1,1 -k 2,2n \
| pigz -c > #{t.name}
EOF
end

##

file 'tmp/byTFE/peaks_nonClass0.bed.gz' => 'out/byTFE/peaks.bed.gz' do |t|
  sh "cp #{t.source} #{t.name}"
end

##

step4a_class2location = ['',
                         "Coding 5'-UTR",
                         "Coding upstream",
                         "Coding CDS",
                         "Coding 3'-UTR",
                         "Noncoding 1st-exon",
                         "Noncoding upstream",
                         "Noncoding other-exon",
                         "Intron",
                         "Unannotated"]

def step4a_peakClass_sources(t)
  cls = /lass(\d+)\.(txt|bed)\.gz$/.match(t).to_a[1].to_i
  return ["#{DEFAULTS['TRANSCRIPT']}.class#{cls}.bed.gz",
          t.sub(/_(nonC|c)lass\d+\.(txt|bed)\.gz$/, "_nonClass#{cls-1}.bed.gz")]
end

rule /peaks_class\d+\.txt\.gz$/ => [->(path){ step4a_peakClass_sources(path) }] do |t|
  cls = /lass(\d+)\.(txt|bed)\.gz$/.match(t.name).to_a[1].to_i
  loc = step4a_class2location[cls]
  
  infp = open("| intersectBed -s -wa -wb -a #{t.source} -b #{t.sources[1]} | gcut -f 1,4,6,9,10")
  tfe2sym2acc = Hash.new
  tfe2id = Hash.new
  while line = infp.gets
    chr, symacc, str, pos, tfe = line.rstrip.split(/\t/)
    sym, acc = symacc.split(/\|/)
    sym = acc if sym == ''
    tfe2sym2acc[tfe] = Hash.new unless tfe2sym2acc.key?(tfe)
    tfe2sym2acc[tfe][sym] = Hash.new unless tfe2sym2acc[tfe].key?(sym)
    tfe2sym2acc[tfe][sym][acc] = ''
    tfe2id[tfe] = "#{chr}\t#{pos}\t#{str}" unless tfe2id.key?(tfe)
  end
  infp.close
  
  outfp = open("| gzip -c > #{t.name}", 'w')
  tfe2sym2acc.each do |tfe, sym2acc|
    syms = Array.new
    symaccs = Array.new
    sym2acc.keys.sort.each do |sym|
      syms.push(sym)
      symaccs.push(sym2acc[sym].keys.sort.join(';'))
    end
    outfp.puts [tfe, tfe2id[tfe], syms.join(';'), symaccs.join('|'), loc].join("\t")
  end
  outfp.close
end

rule /peaks_nonClass[1-8]\.bed\.gz$/ => [->(path){ step4a_peakClass_sources(path) }] do |t|
  cls = /lass(\d+)\.(txt|bed)\.gz$/.match(t.name).to_a[1].to_i
  sh "intersectBed -s -v -b #{t.source} -a #{t.sources[1]} | gzip -c > #{t.name}"
end
                           
##

file 'tmp/byTFE/peaks_class9.txt.gz' =>
     ['out/byTFE/regions.bed.gz', 'tmp/byTFE/peaks_nonClass8.bed.gz'] do |t|
  rid = Hash.new
  infp = open("| gunzip -c #{t.source}")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    rid[cols[3]] = "#{cols[0]}:#{cols[1].to_i+1}-#{cols[2]};#{cols[5]}"
  end
  infp.close
  
  outfp = open("| gzip -c > #{t.name}", 'w')
  infp = open("| gunzip -c #{t.sources[1]}")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    outfp.puts [cols[3], "#{cols[0]}\t#{cols[2]}\t#{cols[5]}", rid[cols[3]], '', 'Unannotated'].join("\t")
  end
  infp.close
  outfp.close
end

##

step4a_annotation_sources = Array.new
1.upto(9) do |i|
  step4a_annotation_sources.push("tmp/byTFE/peaks_class#{i}.txt.gz")
end

file 'out/byTFE/annotation.txt.gz' => step4a_annotation_sources do |t|
  sh "gunzip -c #{t.sources.join(' ')} | gzip -c > #{t.name}"
end
