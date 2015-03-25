##
## Step 2e - base for counting
##

def step2e_sources(path)
  return [path.sub(/^tmp/, 'out/bam').sub(/\.step2e$/, '.bam')]
end

rule /.step2e$/ => [->(path){ step2e_sources(path) }] do |t|
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

#

rule '.step2e_cnt' => '.step2e' do |t|
  sh <<EOF
(gunzip -c #{t.source} | wc -l | gtr -d ' ';\
 gunzip -c #{t.source} | grep ^RNA_SPIKE_ | wc -l | gtr -d ' ')\
| gpaste -s - > #{t.name}
EOF
end

#

task :clean_step2e do
  LIBIDS.each do |libid|
    sh "rm -rf tmp/#{libid}.*.step2e"
    sh "rm -rf tmp/#{libid}.*.step2e_cnt"
  end
end
