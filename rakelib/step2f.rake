##
## Step 2f - Decollapse UMIs
##

def step2f_sources(path)
  return [path.sub(/\.step2f$/, '.step2d'),
          "tmp/#{path.pathmap('%n').pathmap('%n')}.step2a.trace"]
end

rule /\.step2f$/ => [->(path){ step2f_sources(path) }] do |t|
  repacc2accpath = mymkfifo('step2f-')
  pid = spawn "unpigz -c #{t.sources[1]} > #{repacc2accpath}"

  sh <<EOF
unpigz -c #{t.source} \
| gcut -f 2 \
| gsort -S #{80/THREADS}% -t '\t' -k 1,1 \
| gjoin -t '\t' -j 1 - #{repacc2accpath} \
| pigz -c > #{t.name}
EOF

  Process.waitpid(pid)
end

#

rule '.step2f_cnt' => '.step2f' do |t|
  sh "unpigz -c #{t.source} | gwc -l | gtr -d ' ' > #{t.name}"
end
