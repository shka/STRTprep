##
## Step 1b - QV filtering & sort
##

if DEFAULTS['PHYX'].nil?
  def step2b_sources(path)
    libid, runid = parse_librunid(path)
    conf = CONF[libid]
    return conf['FASTQS'][runid]
  end

  rule '.step1b' => [->(path){ step2b_sources(path) }] do |t|
    libid, runid = parse_librunid(t.name)
    conf = CONF[libid]
    len = conf['UMI']+conf['BARCODE']+conf['GAP']+conf['CDNA']
    sh <<EOF
(unpigz -c #{t.source} \
 | ruby -e 'while l=gets; m=gets.rstrip; o=gets; n=gets.rstrip[0, #{len}]; puts "#{libid}\t\#{l.rstrip[1..-1]}:1-#{len}\t\#{n}\t\#{m[0, #{len}]}" if n.index(/[!-#]/).nil?; end' \
 | gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -t '\t' -k 4,4 -k 3,3r \
 | pigz -c > #{t.name}) 2> #{t.name}.log
EOF
  end
else
  rule '.step1b' => '.step1a' do |t|
    libid, runid = parse_librunid(t.name)
    conf = CONF[libid]
    len = conf['UMI']+conf['BARCODE']+conf['GAP']+conf['CDNA']
    sh <<EOF
(samtools view -f 4 -F 256 #{t.source}\
 | gcut -f 1,10,11 \
 | gawk 'BEGIN { FS="\t"; OFS="\t" }; match(substr($3, 1, #{len}), /[!-#]/) == 0 { print "#{libid}",$1":1-#{len}",substr($3, 1, #{len}),substr($2, 1, #{len}) }' \
 | gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -t '\t' -k 4,4 -k 3,3r \
 | pigz -c > #{t.name}) 2> #{t.name}.log
EOF
  end
end
