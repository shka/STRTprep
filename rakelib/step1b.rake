##
## Step 1b - QV filtering & sort
##

rule '.step1b' => '.step1a' do |t|
  libid, runid = parse_librunid(t.name)
  conf = CONF[libid]
  len = conf['UMI']+conf['BARCODE']+conf['GAP']+conf['CDNA']
  sh <<EOF
samtools view -f 4 -F 256 #{t.source}\
| gcut -f 1,10,11\
| gawk 'BEGIN { OFS="\t" }; match(substr($3, 1, #{len}), /[!-#]/) == 0 { print "#{libid}",$1":1-#{len}",substr($3, 1, #{len}),substr($2, 1, #{len}) }'\
| gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -k 4,4 -k 3,3r\
| pigz -c > #{t.name} 2> #{t.name}.log
EOF
end

task :clean_step1b do 
  LIBIDS.each do |libid|
    CONF[libid]['FASTQS'].each_index do |runid|
      rm_rf "tmp/#{libid}.#{runid}.step1b"
    end
  end
end
