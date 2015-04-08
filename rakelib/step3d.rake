##
## Step 3d - merge the counts
##

require 'csv'

step3d_sources = ['src/samples.csv']
LIBWELLIDS.each do |libwellid|
  step3d_sources.push("tmp/byGene/#{libwellid}.step3b")
end

def step3d_job(t, colname)
  libwellid2name = Hash.new
  samples = CSV.table(t.source)
  samples.each do |row|
    libwellid2name["#{row[:library]}.#{row[:well]}"] = row[:name]
  end

  header = [colname]
  LIBWELLIDS.each do |libwellid|
    libwellid2name[libwellid] = 'NA' unless libwellid2name.key?(libwellid)
    header.push("#{libwellid}|#{libwellid2name[libwellid]}")
  end
  
  outfp = open("| gzip -c > #{t.name}", 'w')
  outfp.puts header.join("\t")
  infp = open("| cat "+t.sources[1..-1].join(" | join -t '\t' -j 1 - "))
  while line = infp.gets
    outfp.puts line
  end
  infp.close
  outfp.close
end

file 'out/byGene/reads_all.txt.gz' => step3d_sources do |t|
  step3d_job(t, 'Gene')
end

##

task :clean_step3d do
  sh 'rm out/byGene/reads_all.txt.gz'
end
