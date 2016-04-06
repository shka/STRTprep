##
## Step 5 - files for exports
##

##

file 'out/web/regions_byGene.bed.gz' => 'out/byGene/regions.bed.gz' do |t|
  mkdir_p t.name.pathmap('%d')
  pre = t.name.pathmap('%X')
  sh "echo 'track name=regions_byGene description=\"Regions for gene-based quantification\" visibility=pack colorByStrand=\"133,153,0 203,75,22\"' > #{pre}"
  sh "unpigz -c #{t.source} | grep ^chr >> #{pre}"
  sh "rm -f #{t.name}"
  sh "pigz #{pre}"
end

##

rule /_transcripts\.gtf\.gz$/ => [->(path){ path.sub(/^out\/web\//, 'tmp/byTFE/').sub(/_transcripts\.gtf\.gz$/, '.step4a/transcripts.gtf') }] do |t|
  mkdir_p t.name.pathmap('%d')
  pre = t.name.pathmap('%X')
  cls = /\/([^\/]+)_transcripts\.gtf\.gz$/.match(t.name).to_a[1]
  sh "echo 'track name=#{cls} descriptsion=\"Assembled transcripts by #{cls} samples\" visibility=pack' > #{pre}"
  sh "grep ^chr #{t.source} >> #{pre}"
  sh "rm -f #{t.name}"
  sh "pigz #{pre}"
end
                         
##

file 'out/web/regions_byTFE.bed.gz' => ['out/byTFE/annotation.txt.gz', 'out/byTFE/regions.bed.gz'] do |t|
  peaks = Hash.new
  infp = open("| unpigz -c #{t.source}")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    peaks[cols[0]] = cols[2].to_i
  end
  infp.close

  mkdir_p t.name.pathmap('%d')
  pre = t.name.pathmap('%X')
  sh "echo 'track name=regions_byTFE description=\"Regions for TFE-based quantification, and the peaks\" visibility=pack colorByStrand=\"133,153,0 203,75,22\"' > #{pre}"
  outfp = open(pre, 'a')
  infp = open("| unpigz -c #{t.sources[1]} | grep ^chr")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    if peaks.key?(cols[3])
      peak = peaks[cols[3]]
      outfp.puts (cols + [peak-1, peak]).join("\t")
    end
  end
  infp.close
  outfp.close
  sh "rm -f #{t.name}"
  sh "pigz #{pre}"
end
##

require 'csv'

step5_targets = ['out/web/regions_byTFE.bed.gz']
begin
  tmp = Hash.new
  samples = CSV.table('src/samples.csv')
  samples.each do |row|
    cls = row[:classtfe]
    tmp[cls] = '' if cls != 'NA' && cls != ''
  end
  tmp.each_key do |cls|
    step5_targets.push("out/web/#{cls}_transcripts.gtf.gz")
  end
end

task :web => step5_targets
