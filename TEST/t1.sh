../ONEview -b t1.foo > ZZ_1.1foo
../ONEview -s ZZ_1.1foo > ZZ_1.schema
../ONEview -h ZZ_1.1foo > ZZ_1.body
../ONEview -S ZZ_1.schema -t foo ZZ_1.body
