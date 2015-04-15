##
## Step 3g - fluctuated genes
##

file 'out/byGene/fluctuation.txt.gz' => 'out/byGene/nreads.RData' do |t|
  sh "R --vanilla --quiet --args #{DEFAULTS['FLUCTUATION']} #{t.name.pathmap('%d')} < bin/_step3g_fluctuation.R > #{t.name}.log 2>&1"
end

task :clean_step3g do
  rm 'out/byGene/fluctuation.RData'
  rm 'out/byGene/fluctuation.txt.gz'
end

