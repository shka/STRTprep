##
## Step 4d - merge the counts
##

require 'csv'

step4d_sources = ['src/samples.csv']
LIBWELLIDS.each do |libwellid|
  step4d_sources.push("tmp/byTFE/#{libwellid}.step4b")
end

file 'out/byTFE/reads_all.txt.gz' => step4d_sources do |t|
  step3d_job(t, 'TFE')
end

##

task :clean_step4d do
  sh 'rm out/byTFE/reads_all.txt.gz'
end
