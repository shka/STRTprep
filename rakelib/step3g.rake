##
## Step 3g - fluctuated genes
##

file 'out/cg/fluctuation.RData' => 'out/cg/nreads.RData' do |t|
  sh "R --vanilla --quiet < bin/_step3g_fluctuation.R > #{t.name}.log 2>&1"
end

file 'out/cg/fluctuation.txt.gz' => 'out/cg/fluctuation.RData'

task :clean_step3g do
  rm 'out/cg/fluctuation.RData'
  rm 'out/cg/fluctuation.txt.gz'
end
