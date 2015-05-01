##
## Step 4h - differential expression analysis
##

step4h_sources = ['out/byTFE/reads.txt.gz',
                  'out/byTFE/nreads.txt.gz',
                  'out/byTFE/fluctuation.txt.gz',
                  'out/byTFE/regions.bed.gz',
                  'out/byTFE/annotation.txt.gz']
begin
  infp = open('src/samples.csv', 'rt')
  colnames = infp.gets.rstrip.split(',')
  classes = colnames.select { |colname| /^CLASS\.\d+$/.match(colname) }
  classes.each do |cls|
    idx = /^CLASS\.(\d+)$/.match(cls).to_a[1]
    step4h_sources.push("out/byTFE/diffexp#{idx}.txt.gz")
  end
  infp.close
end

file 'out/byTFE/diffexp.csv' => step4h_sources do |t|
  step3h_job(t)
end
