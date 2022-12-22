mkdir out

/usr/bin/time -v python3 conv_study.py --example 2 --n_ref 5 --order 2 --struc_mesh 1 -def 0 --method trefftz
/usr/bin/time -v python3 conv_study.py --example 2 --n_ref 5 --order 2 --struc_mesh 1 -def 0 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 2 --n_ref 5 --order 2 --struc_mesh 1 -def 0 --method l2

/usr/bin/time -v python3 conv_study.py --example 2 --n_ref 5 --order 3 --struc_mesh 1 -def 0 --method trefftz
/usr/bin/time -v python3 conv_study.py --example 2 --n_ref 5 --order 3 --struc_mesh 1 -def 0 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 2 --n_ref 5 --order 3 --struc_mesh 1 -def 0 --method l2

/usr/bin/time -v python3 conv_study.py --example 2 --n_ref 4 --order 4 --struc_mesh 1 -def 0 --method trefftz
/usr/bin/time -v python3 conv_study.py --example 2 --n_ref 4 --order 4 --struc_mesh 1 -def 0 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 2 --n_ref 4 --order 4 --struc_mesh 1 -def 0 --method l2