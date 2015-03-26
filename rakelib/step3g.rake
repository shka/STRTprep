##
## Step 3g - fluctuated genes
##

rule 'out/cg/fluctuation.RData' => 'out/cg/nreads.RData' do |t|
  sh "R --vanilla --quiet < bin/_step3g_fluctuation.R > #{t.name}.log 2>&1"
end
