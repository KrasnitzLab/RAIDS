
############################################################
## How to create the pedSyn.RDS file
## This file contains a small synthetic pedigree data.frame
############################################################

synthetic <- data.frame(data.id=c('1.ex1.HG00246.1', '1.ex1.HG00138.1',
    '1.ex1.HG00325.1', '1.ex1.HG00330.1', '1.ex1.HG00611.1',
    '1.ex1.HG00406.1', '1.ex1.HG01173.1', '1.ex1.HG01403.1',
    '1.ex1.HG02165.1', '1.ex1.HG02390.1', '1.ex1.HG01112.1',
    '1.ex1.HG01142.1', '1.ex1.HG01615.1', '1.ex1.HG01606.1',
    '1.ex1.HG01968.1', '1.ex1.HG02299.1', '1.ex1.HG02658.1',
    '1.ex1.HG03238.1', '1.ex1.HG01850.1', '1.ex1.HG02016.1',
    '1.ex1.HG02013.1', '1.ex1.HG02339.1', '1.ex1.HG02465.1',
    '1.ex1.HG02884.1', '1.ex1.HG02974.1', '1.ex1.HG03300.1',
    '1.ex1.HG03814.1', '1.ex1.HG04153.1', '1.ex1.HG03445.1',
    '1.ex1.HG03557.1', '1.ex1.HG03689.1', '1.ex1.HG03848.1',
    '1.ex1.HG03789.1', '1.ex1.HG04239.1', '1.ex1.NA12751.1',
    '1.ex1.NA12717.1', '1.ex1.NA19107.1', '1.ex1.NA19143.1',
    '1.ex1.NA18548.1', '1.ex1.NA18616.1', '1.ex1.NA19075.1',
    '1.ex1.NA19010.1', '1.ex1.NA19475.1', '1.ex1.NA19466.1',
    '1.ex1.NA19712.1', '1.ex1.NA20126.1', '1.ex1.NA19731.1',
    '1.ex1.NA19794.1', '1.ex1.NA20528.1', '1.ex1.NA20757.1',
    '1.ex1.NA20908.1', '1.ex1.NA20847.1'),
    case.id=c('HG00246', 'HG00138', 'HG00325', 'HG00330', 'HG00611', 'HG00406',
    'HG01173', 'HG01403', 'HG02165', 'HG02390', 'HG01112', 'HG01142',
    'HG01615', 'HG01606', 'HG01968', 'HG02299', 'HG02658', 'HG03238',
    'HG01850', 'HG02016', 'HG02013', 'HG02339', 'HG02465', 'HG02884',
    'HG02974', 'HG03300', 'HG03814', 'HG04153', 'HG03445', 'HG03557',
    'HG03689', 'HG03848', 'HG03789', 'HG04239', 'NA12751', 'NA12717',
    'NA19107', 'NA19143', 'NA18548', 'NA18616', 'NA19075', 'NA19010',
    'NA19475', 'NA19466', 'NA19712', 'NA20126', 'NA19731', 'NA19794',
    'NA20528', 'NA20757', 'NA20908', 'NA20847'),
    sample.type=rep('Synthetic', 52),
    diagnosis=rep('Cancer', 52),
    source=rep('Synthetic', 52),
    study.id=rep('MYDATA.Synthetic', 52),
    superPop=c('EUR', 'EUR', 'EUR', 'EUR', 'EAS', 'EAS', 'AMR', 'AMR', 'EAS',
    'EAS', 'AMR', 'AMR', 'EUR', 'EUR', 'AMR', 'AMR', 'SAS', 'SAS', 'EAS',
    'EAS', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'AFR', 'SAS', 'SAS', 'AFR',
    'AFR', 'SAS', 'SAS', 'SAS', 'SAS', 'EUR', 'EUR', 'AFR', 'AFR', 'EAS',
    'EAS', 'EAS', 'EAS', 'AFR', 'AFR', 'AFR', 'AFR', 'AMR', 'AMR', 'EUR',
    'EUR', 'SAS', 'SAS'), stringsAsFactors=FALSE)

rownames(synthetic) <- synthetic$data.id

saveRDS(synthetic, file="pedSyn.RDS")
