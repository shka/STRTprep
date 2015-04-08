##
## Step 4g - fluctuated genes
##

file 'out/byTFE/fluctuation.txt.gz' => 'out/byTFE/nreads.RData' do |t|
  sh "R --vanilla --quiet --args #{DEFAULTS['FLUCTUATION']} #{t.name.pathmap('%d')} < bin/_step3g_fluctuation.R > #{t.name}.log 2>&1"
end
