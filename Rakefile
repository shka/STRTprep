load 'rake/miscellaneous.rake'
load 'rake/configuration.rake'
load 'rake/removePhyX.rake'
load 'rake/demultiplex.rake'
load 'rake/alignment.rake'
load 'rake/hub.rake'
load 'rake/expression.rake'

task :default => :expression

task 'clean' => ['clean_removePhyX', 'clean_demultiplex', 'clean_alignment', 'clean_hub', 'clean_expression']
