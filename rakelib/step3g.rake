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
      `md5sum #{t.source} #{tmp} | gcut -d ' ' -f 1 | guniq | gwc -l`.to_i != 1)
    sh "bin/_step3g_fluctuation.R #{DEFAULTS['FLUCTUATION']} #{dir} > #{t.name}.log 2>&1"
    sh "cp -p #{t.source} #{tmp}"
  else
    puts "... skipped #{t.name}, since the qualified samples were identical with the previous run"
    sh "touch #{t.name}"
  end
end

##

def step3g_fluctuationCSV_source(path)
  return path.sub(/fluctuation\.csv/, 'diffexp.csv')
end

rule /\/fluctuation\.csv$/ => [->(p){ step3g_fluctuationCSV_source(p) }] do |t|
  outfp = open(t.name, 'w')
  infp = open(t.source)
  headers = infp.gets.rstrip.split(/,/)
  targets = Array.new
  0.upto(headers.length-1) do |i|
    targets.push(i) if /^(Score|pvalue|qvalue|fluctuation)\.\d+$/.match(headers[i]).nil?
  end
  headers[headers.index('fluctuation.global')] = 'fluctuation'
  outfp.puts (targets.map { |i| headers[i] }).join(',')
  while line = infp.gets
    cols = line.rstrip.split(/,/)
    outfp.puts (targets.map { |i| cols[i] }).join(',')
  end
  infp.close
  outfp.close
end

##

task :step3g => 'out/byGene/fluctuation.xls'
