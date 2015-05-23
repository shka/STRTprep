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
    mkdir_p t.name.pathmap('%d')
    sh <<EOF
(unpigz -c #{t.source} \
 | bin/_step1b_fastq2oneLine \
 | bin/_step1b_trimWithQCFilter #{libid} #{len} #{t.name}.stat \
 | gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -t '\t' -k 4,4 -k 3,3r \
 | pigz -c > #{t.name}) 2> #{t.name}.log
EOF
    sh "touch #{t.name}.stat"
  end
else
  rule '.step1b' => '.step1a' do |t|
    libid, runid = parse_librunid(t.name)
    conf = CONF[libid]
    len = conf['UMI']+conf['BARCODE']+conf['GAP']+conf['CDNA']
    sh <<EOF
(samtools view -f 4 -F 256 #{t.source}\
 | gcut -f 1,10,11 \
 | bin/_step1b_trimWithQCFilter #{libid} #{len} #{t.name}.stat \
 | gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -t '\t' -k 4,4 -k 3,3r \
 | pigz -c > #{t.name}) 2> #{t.name}.log
EOF
    sh "touch #{t.name}.stat"
  end
end

rule /\.step1b\.stat$/ => [->(p){ p.pathmap("%X") }]

def step1b_byLibraries_sources(path)
  tmp = Array.new
  libid = path.pathmap("%X").pathmap("%X").pathmap("%f").sub('fig.', '')
  CONF[libid]['FASTQS'].each_index do |runid|
    tmp.push("tmp/#{libid}.#{runid}.step1b.stat")
  end
  return tmp
end

rule /\.qualifiedReads\.pdf$/ => [->(p){ step1b_byLibraries_sources(p) }] do |t|
  libid = t.name.pathmap("%X").pathmap("%X").pathmap("%f").sub('fig.', '')
  conf = CONF[libid]
  len = conf['UMI']+conf['BARCODE']+conf['GAP']+conf['CDNA']
  mkdir_p 'out'
  sh "bin/_step1b_qualifiedReads.R #{libid} #{len} #{t.sources.join(' ')} > #{t.name}.log 2>&1"
end
