##
## Step 3b - base for counting
##

def step3b_sources(path)
  return [path.sub(/^tmp/, 'out/bam').sub(/\.step3b$/, '.bam')]
end

rule /^tmp\/[^\/]+\.step3b$/ => [->(path){ step3b_sources(path) }] do |t|
  outfp = open("| gzip -c > #{t.name}", 'w')
  infp = open("| bamToBed -i #{t.source}")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    acc = "#{cols[3]}:5end"
    if cols[5] == '+'
      outfp.puts [cols[0], cols[1], cols[1].to_i+1, acc, cols[4], cols[5]].join("\t")
    else
      outfp.puts [cols[0], cols[2].to_i-1, cols[2], acc, cols[4], cols[5]].join("\t")
    end
  end
  infp.close
  outfp.close
end

task :clean_step3b do
  LIBIDS.each do |libid|
    sh "rm -rf tmp/#{libid}.*.step3b"
  end
end
