##
## Step 3g - fluctuated genes
##

def step3g_fluctuation_source(path)
  return "#{path.pathmap('%d')}/nreads.RData"
end

rule /\/fluctuation\.txt\.gz$/ => [->(p){ step3g_fluctuation_source(p) }] do |t|
  dir = t.name.pathmap('%d')
  tmp = "#{dir.sub(/^out/, 'tmp')}/nreads.RData"
  if (!File.exist?(t.name) ||
      !File.exist?(tmp) ||
      `gmd5sum #{t.source} #{tmp} | gcut -d ' ' -f 1 | guniq | gwc -l`.to_i != 1)
    sh "bin/_step3g_fluctuation.R #{dir} > #{t.name}.log 2>&1"
    sh "cp -p #{t.source} #{tmp}"
  else
    puts "... skipped #{t.name}, since the qualified samples were identical with the previous run"
    sh "touch #{t.name}"
  end
end
