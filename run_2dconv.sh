mkdir out

/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 2 -def 0 --method trefftz
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 2 -def 0 --method trefftz --stabil patch
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 2 -def 0 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 2 -def 0 --method l2
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 2 -def 0 --method l2 --stabil patch

/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 3 -def 0 --method trefftz
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 3 -def 0 --method trefftz --stabil patch 
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 3 -def 0 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 3 -def 0 --method l2
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 7 --order 3 -def 0 --method l2 --stabil patch

/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 5 --order 4 -def 0 --method trefftz
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 5 --order 4 -def 0 --method trefftz --stabil patch
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 5 --order 4 -def 0 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 5 --order 4 -def 0 --method l2
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 5 --order 4 -def 0 --method l2 --stabil patch

/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 4 --order 5 -def 0 --method trefftz
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 4 --order 5 -def 0 --method trefftz --stabil patch
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 4 --order 5 -def 0 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 4 --order 5 -def 0 --method l2
/usr/bin/time -v python3 conv_study.py --example 1 --n_ref 4 --order 5 -def 0 --method l2 --stabil patch


/usr/bin/time -v python3 conv_study.py --example 3 --n_ref 7 --order 2 -def 1 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 3 --n_ref 7 --order 2 -def 1 --method l2

/usr/bin/time -v python3 conv_study.py --example 3 --n_ref 7 --order 3 -def 1 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 3 --n_ref 7 --order 3 -def 1 --method l2

/usr/bin/time -v python3 conv_study.py --example 3 --n_ref 7 --order 4 -def 1 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 3 --n_ref 7 --order 4 -def 1 --method l2

/usr/bin/time -v python3 conv_study.py --example 3 --n_ref 6 --order 5 -def 1 --method emb_trefftz
/usr/bin/time -v python3 conv_study.py --example 3 --n_ref 6 --order 5 -def 1 --method l2