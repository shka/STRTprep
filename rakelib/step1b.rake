##
## Step 1b - QV filtering & sort
##

qb = !PREPROCESS.key?('QUALITYBASE') || PREPROCESS['QUALITYBASE'] == 33 ? 33 : PREPROCESS['QUALITYBASE']

if !PREPROCESS.key?('PHYX') || PREPROCESS['PHYX'].nil?

  if PREPROCESS.key?('MULTIPLEXTYPE') && PREPROCESS['MULTIPLEXTYPE'] == 'C1'

    def step1b_sources(path)
      libid, runid = parse_librunid(path)
      return [LIBRARIES[libid]['FASTQS'][runid],
              LIBRARIES[libid]['FASTQS2'][runid]]
    end

    rule '.step1b' => [->(path){ step1b_sources(path) }] do |t|
      libid, runid = parse_librunid(t.name)
      len = PREPROCESS['UMI']+PREPROCESS['BARCODE']+PREPROCESS['GAP']+PREPROCESS['CDNA']
      mkdir_p t.name.pathmap('%d')
      read1fifo = mymkfifo('step1b-read1-')
      pid1 = spawn "unpigz -c #{t.sources[0]} > #{read1fifo}"
      read2fifo = mymkfifo('step1b-read2-')
      pid2 = spawn "unpigz -c #{t.sources[1]} > #{read2fifo}"
      sh <<EOF
(gpaste #{read1fifo} #{read2fifo} \
 | bin/_step1b_fastqs2oneLine #{PREPROCESS['UMI']} \
 | bin/_step1b_trimWithQCFilter #{libid} #{len} #{t.name}.stat #{qb} \
 | gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -t '\t' -k 4,4 -k 3,3r \
 | pigz -c > #{t.name}) 2> #{t.name}.log
EOF
      sh "touch #{t.name}.stat"
      Process.waitpid(pid1)
      Process.waitpid(pid2)
    end
    
  else
    
    def step1b_sources(path)
      libid, runid = parse_librunid(path)
      return LIBRARIES[libid]['FASTQS'][runid]
    end
    
    rule '.step1b' => [->(path){ step1b_sources(path) }] do |t|
      libid, runid = parse_librunid(t.name)
      len = PREPROCESS['UMI']+PREPROCESS['BARCODE']+PREPROCESS['GAP']+PREPROCESS['CDNA']
      mkdir_p t.name.pathmap('%d')
      sh <<EOF
(unpigz -c #{t.source} \
 | bin/_step1b_fastq2oneLine \
 | bin/_step1b_trimWithQCFilter #{libid} #{len} #{t.name}.stat #{qb} \
 | gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -t '\t' -k 4,4 -k 3,3r \
 | pigz -c > #{t.name}) 2> #{t.name}.log
EOF
      sh "touch #{t.name}.stat"
    end
    
  end
  
else
  
  rule '.step1b' => '.step1a' do |t|
    libid, runid = parse_librunid(t.name)
    len = PREPROCESS['UMI']+PREPROCESS['BARCODE']+PREPROCESS['GAP']+PREPROCESS['CDNA']
    sh <<EOF
(samtools view -f 4 -F 256 #{t.source}\
 | gcut -f 1,10,11 \
 | bin/_step1b_trimWithQCFilter #{libid} #{len} #{t.name}.stat #{qb} \
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
  LIBRARIES[libid]['FASTQS'].each_index do |runid|
    tmp.push("tmp/#{libid}.#{runid}.step1b.stat")
  end
  return tmp
end

rule /\.qualifiedReads\.pdf$/ => [->(p){ step1b_byLibraries_sources(p) }] do |t|
  libid = t.name.pathmap("%X").pathmap("%X").pathmap("%f").sub('fig.', '')
  len = PREPROCESS['UMI']+PREPROCESS['BARCODE']+PREPROCESS['GAP']+PREPROCESS['CDNA']
  mkdir_p 'out'
  sh "bin/_step1b_qualifiedReads.R #{libid} #{len} #{t.sources.join(' ')} > #{t.name}.log 2>&1"
end
