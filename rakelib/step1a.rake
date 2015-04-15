##
## Step 1a - exclusion of PhyX reads
##

unless DEFAULTS['PHYX'].nil?
  def step1a_sources(path)
    libid, runid = parse_librunid(path)
    conf = CONF[libid]
    return [ conf['FASTQS'][runid],
             File.expand_path(DEFAULTS['PHYX']+'.1.ebwt') ]
  end
  
  rule '.step1a' => [->(path){ step1a_sources(path) }] do |t|
    mkdir_p 'tmp'
    sh <<EOF
(gunzip -c #{t.source}\
| bowtie -S -p #{PROCS} -t #{t.sources[1].pathmap('%X').pathmap('%X')} -\
| samtools view -S -b -o #{t.name} -) > #{t.name}.log 2>&1
EOF
  end
end
