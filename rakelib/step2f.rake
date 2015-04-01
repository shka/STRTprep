##
## Step 2f - Decollapse UMIs
##

def step2f_sources(path)
  return [path.sub(/\.step2f$/, '.step2d'), 'tmp/step2a.trace']
end

rule /\.step2f$/ => [->(path){ step2f_sources(path) }] do |t|
  repacc2accpath = mymkfifo('step2f-')
  pid = spawn "gunzip -c #{t.sources[1]} > #{repacc2accpath}"

  sh <<EOF
gunzip -c #{t.source} \
| gcut -f 2 \
| gsort -S 3G -k 1,1\
| gjoin -t '\t' -j 1 - #{repacc2accpath}\
| gzip -c > #{t.name}
EOF

  Process.waitpid(pid)
end

#

rule '.step2f_cnt' => '.step2f' do |t|
  sh "gunzip -c #{t.source} | wc -l | gtr -d ' ' > #{t.name}"
end

#

task :clean_step2f do |t|
  LIBIDS.each do |libid|
    sh "rm tmp/#{libid}.*.step2f"
    sh "rm tmp/#{libid}.*.step2f_cnt"
  end
end
