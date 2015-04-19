##
## Step 2e - base for counting
##

def step2e_sources(path)
  return [path.sub(/^tmp/, 'out/bam').sub(/\.step2e$/, '.bam')]
end

rule /.step2e$/ => [->(path){ step2e_sources(path) }] do |t|
  sh <<EOF
bamToBed -i #{t.source} \
| ruby -anle 'puts "\#{$F[0]}\t\#{$F[5]=="+"?$F[1]:$F[2].to_i-1}\t\#{$F[5]=="+"?$F[1].to_i+1:$F[2]}\t\#{$F[3]}:5end\t\#{$F[4]}\t\#{$F[5]}"' \
| pigz -c > #{t.name}
EOF
end

#

rule '.step2e_cnt' => '.step2e' do |t|
  sh <<EOF
(unpigz -c #{t.source} | wc -l | gtr -d ' '; \
 unpigz -c #{t.source} | grep ^RNA_SPIKE_ | wc -l | gtr -d ' ')\
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
