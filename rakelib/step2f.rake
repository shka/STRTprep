##
## Step 2f - Decollapse UMIs
##

def step2f_sources(path)
  return [path.sub(/\.step2f$/, '.step2d'), 'tmp/step2a.trace']
end

rule /.step2f$/ => [->(path){ step2f_sources(path) }] do |t|
  repacc2accpath = mymkfifo('step2f-')
  pid1 = spawn "gunzip -c #{t.sources[1]} > #{repacc2accpath}"

  libwellid = t.name.pathmap("%n")
  repaccpath = mymkfifo('step2f-')
  pid2 = spawn <<EOF
gunzip -c #{t.source} \
| gcut -f 2 \
| gsort -S #{1500*PROCS}M -k 1,1 > #{repaccpath}
EOF

  sh <<EOF
gjoin -t '\t' -j 1 #{repaccpath} #{repacc2accpath}\
| gsort -S #{1500*PROCS}M -k 1,1\
| gzip -c > #{t.name}
EOF

  Process.waitpid(pid1)
  Process.waitpid(pid2)
end

#

rule '.step2f_cnt' => '.step2f' do |t|
  sh "gunzip -c #{t.source} | wc -l | gtr -d ' ' > #{t.name}"
end
