V33 :0x4 initialize
18 mod_initialize.f90 S624 0
11/27/2019  13:06:51
use units public 0 indirect
use setup_variables public 0 indirect
use route_control public 0 indirect
use global_variables public 0 direct
enduse
S 624 24 0 0 0 6 1 0 5015 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 initialize
S 1415 23 5 0 0 0 1416 624 11765 0 0 A 0 0 0 0 B 0 118 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 read_input
S 1416 14 5 0 0 0 1 1415 11765 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 23 0 0 0 0 0 0 0 0 0 0 0 0 0 24 0 624 0 0 0 0 read_input
F 1416 0
S 1417 23 5 0 0 0 1418 624 11776 0 0 A 0 0 0 0 B 0 181 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 set_default
S 1418 14 5 0 0 0 1 1417 11776 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 122 0 624 0 0 0 0 set_default
F 1418 0
S 1419 23 5 0 0 0 1420 624 11788 0 0 A 0 0 0 0 B 0 231 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 read_tdcihead
S 1420 14 5 0 0 0 1 1419 11788 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 25 0 0 0 0 0 0 0 0 0 0 0 0 0 185 0 624 0 0 0 0 read_tdcihead
F 1420 0
S 1421 23 5 0 0 0 1422 624 11802 0 0 A 0 0 0 0 B 0 501 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 set_variables
S 1422 14 5 0 0 0 1 1421 11802 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 235 0 624 0 0 0 0 set_variables
F 1422 0
S 1423 23 5 0 0 0 1425 624 11816 0 0 A 0 0 0 0 B 0 676 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 allocate_main
S 1424 1 3 1 0 28 1 1423 11830 4 43000 A 0 0 0 0 B 0 676 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 option
S 1425 14 5 0 0 0 1 1423 11816 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 27 1 0 0 0 0 0 0 0 0 0 0 0 0 505 0 624 0 0 0 0 allocate_main
F 1425 1 1424
S 1426 23 5 0 0 0 1428 624 11837 0 0 A 0 0 0 0 B 0 850 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 deallocate_main
S 1427 1 3 1 0 28 1 1426 11830 4 43000 A 0 0 0 0 B 0 850 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 option
S 1428 14 5 0 0 0 1 1426 11837 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 29 1 0 0 0 0 0 0 0 0 0 0 0 0 680 0 624 0 0 0 0 deallocate_main
F 1428 1 1427
S 1429 23 5 0 0 0 1430 624 11853 0 0 A 0 0 0 0 B 0 933 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 read_hamdata
S 1430 14 5 0 0 0 1 1429 11853 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 31 0 0 0 0 0 0 0 0 0 0 0 0 0 854 0 624 0 0 0 0 read_hamdata
F 1430 0
S 1431 23 5 0 0 0 1433 624 11866 0 0 A 0 0 0 0 B 0 1027 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 write_input
S 1432 1 3 1 0 7 1 1431 11878 4 3000 A 0 0 0 0 B 0 1027 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 myoption
S 1433 14 5 0 0 0 1 1431 11866 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 32 1 0 0 0 0 0 0 0 0 0 0 0 0 937 0 624 0 0 0 0 write_input
F 1433 1 1432
S 1434 23 5 0 0 0 1437 624 11887 0 0 A 0 0 0 0 B 0 1082 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 writeme_modreadin
S 1435 1 3 1 0 28 1 1434 11905 4 43000 A 0 0 0 0 B 0 1082 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 myroutine
S 1436 1 3 1 0 28 1 1434 11830 4 43000 A 0 0 0 0 B 0 1082 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 option
S 1437 14 5 0 0 0 1 1434 11887 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 34 2 0 0 0 0 0 0 0 0 0 0 0 0 1031 0 624 0 0 0 0 writeme_modreadin
F 1437 2 1435 1436
S 1438 23 5 0 0 0 1441 624 11915 0 0 A 0 0 0 0 B 0 1154 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 errors_modreadin
S 1439 1 3 1 0 28 1 1438 11905 4 43000 A 0 0 0 0 B 0 1154 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 myroutine
S 1440 1 3 1 0 28 1 1438 11830 4 43000 A 0 0 0 0 B 0 1154 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 option
S 1441 14 5 0 0 0 1 1438 11915 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 37 2 0 0 0 0 0 0 0 0 0 0 0 0 1086 0 624 0 0 0 0 errors_modreadin
F 1441 2 1439 1440
Z
Z