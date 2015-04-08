##
## Step 4f - Reads and the normalized reads in the qualified samples
##


file 'out/byTFE/reads.txt.gz' =>
     ['out/byGene/samples.csv', 'out/byTFE/reads_all.txt.gz'] do |t|
  step3f_job(t, 'TFE')
end
