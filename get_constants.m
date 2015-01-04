function constants = get_constants()

    constants.celldata_version = 1.86;

    constants.f = @z2ipldepth;
    constants.str = @sa2strength;
    
    constants.res = [16.5 16.5 23]*2*1.15;
    constants.strat_x = -20:120;
    constants.point_dir = '/data/home/greenem/data/stratification/surface_points_trans/';
    constants.raw_point_dir = '/data/home/greenem/data/stratification/point_data/';
    constants.surface_point_dir = '/data/home/greenem/data/stratification/surface_points/';
    
    
    constants.raw_conn_loc = '/data/home/greenem/data/stratification/conns.mat';
    constants.conn_loc = '/data/home/greenem/data/stratification/trans_conns.mat';
    constants.trans_loc = '/data/home/greenem/data/stratification/T.mat';
    
    constants.soma_loc_threshold = 20;
    constants.soma_dir = '/data/home/greenem/data/stratification/soma_points/';
    constants.dist_bin = 10;
    constants.angle_step = 2*pi/100;
    constants.strat_dir = '/data/home/greenem/data/stratification/cell_strats/';
    constants.axial_x_min = -20;
    constants.celldata_dir = '/data/home/greenem/data/stratification/cell_data/';
    
    constants.sac_bins = 10:20:20*8.5;
    constants.j_bins = 15:30:30*6;
    constants.minij_bins = 7.5:15:15*7.5;
%     constants.sac_bins = 5:10:135;
    
    constants.min_xy = [.05 0] * 10^5;
    constants.max_xy = [3.45 3.92] * 10^5;

    constants.hull_dir = '/data/home/greenem/data/stratification/hull_intersects/';
    
    constants.cell_dsmp_fact = [round(abs(constants.f(1000)-constants.f(0))) 1000 1000 ] .* [1 1 1];
    
%     constants.colormap = [ 0 137 68; 234 185 31; 37 170 226; 11 70 121; 146 39 143; 196 29 39; 0 0 0]/255;
    constants.colormap = [0 0 4; 0 2 0; 4 0 0; 0 3 3; 3 0 3; 3 3 0; 1 1 1; 0 0 0; 2 3 2]/4;
    
    constants.type.oodsgc = [17161       17080       90001       90002];
    
    constants.type.ganglion = [10005 10007 10010 10013 10015 10017 10018 15018 15066 17009 17011 17012 ...
        17013 17021 17022 17024 17027 17028 17034 17035 17037 17038 17040 17050 17051 17053 17055 17057 17059 17060 ...
        17061 17062 17064 17069 17071 17075 17076 17077 17078 17079 17080 17081 17082 17083 17084 17090 17092 17093 ...
        17095 17097 17098 17105 17107 17109 17110 17111 17114 17118 17121 17126 17127 17130 17132 17134 17135 17138 ...
        17140 17144 17146 17151 17152 17159 17160 17161 17167 17168 17176 17177 17181 17182 17188 17190 17192 17200 ...
        17205 17212 17216 17236 17238 17247 30001 30002 30003 50001 50002 50004 90001 90002];
    
    constants.type.t1 = unique([60008 60019 60026 60027 60032 60052 60055 60079 60099 60105 60109 60110 60111 60114 60118 60129 60132 60139 60147 60150 60158 60161 60162 60164 60170 60177 60184 60188 60189 60142 60194 60195 60196 60203 60204 60212 60213 60216 60218  ...
        60170 60196 60027 60078 60189 60008 60114 60118 60139 60164 60218 60026 60032 60052 60055 60105 60109 60132 60150 60158 60177 60184 60203 60019 60079 60099 60147 60162 60142 60194 60216 60187 60161 60195 60204 60213 60188]); % 60248 60272 60230 60277 60227 60279]); % 60297 60303 60307]);
    constants.type.t2 = unique([60001 60002 60003 60004 60006 60009 60010 60011 60013 60021 60022 60023 60025 60029 60037 60038 60039 60040 60041 60042 60043 60097 60101 60102 60103 60104 60112 60120 60124 60130 60133 60135 60138 60140 60141 60149 60157 60159 60167 60169 60182 60012 60046 60198 60208 60209 60210 60214 60217 60219 ... %   60224 60226 ...
        60039 60140 60182 60012 60210  60013 60104 60135 60149 60214 60221 60003 60037 60038 60040 60157 60223 60001 60002 60022 60023 60042 60097 60112 60120 60130 60133 60138 60198 60208 60209  60004 60010 60011 60025 60029 60041 60101 60102 60141 60159 60167 60219 60080 60006 60021 60043 60103 60124 60009 60046 60044 60224 60221 60223 60224 60226]); %60226  60235 60236 60241 60242 60245 60254 60257 60261 60263 60266 60267 60268 60269 60270 60274 60285 60286]); %60291 60292 60296 60299 60302 60308 60316 60319 60323 60325 60326 60327]);
    constants.type.old_t3a = [60024 60028 60030 60048 60049 60050 60056 60066 60068 60072 60075 60076 60085 60092 60093 60096 60100 60108 60115 60123 60134 60136 60137 60146 60160 60166 60172 60181 60183 60197];
    constants.type.old_t3b = [60034 60045 60047 60053 60057 60058 60063 60069 60073 60074 60077 60081 60082 60084 60090 60091 60094 60106 60121 60122 60125 60127 60148 60143 60151 60153 60154 60165 60168 60173 60175 60179 60185 60190 60191 60200 60201 60202 60206 60211 60215 60222 60228];
    constants.type.old_t4 = [60059 60067 60083 60086 60087 60088 60089 60095 60098 60107 60116 60117 60119 60131 60144 60145 60156 60163 60171 60174 60176 60178 60186 60199 60220];    
        constants.type.t4 = [89 83 119 163 95 186 160 191 171 80 144 222 117 202 58 67 125 116 168 47 121 183 86 206 175 91 197 100 143 127 185 199 220 174] + 60000;
    constants.type.t3b = [96 151 131 94 90 69 215 45 24 63 228 154 122 81 87 173 201 84 57 106 73 178 156 107 211 165 77 153 82 53 190 148 74 98 179] + 60000;
    constants.type.t3a = [93 92 145 34 108 59 181 137 123 176 115 28 146 136 66 76 172 48 50 166 88 72 68 134 85 49 56 75 30] + 60000;    
    constants.type.A2 = [70012 70039 70044 70049 70051 70052 70053 70054 70056 70057 70058 70059 70060 70061 70062 70063 70064 70065 70072 70074 70078];
    
    
    constants.type.on_bc = [60016 60017 60018 60020 60031 60034 60035 60036 60054 60061 60064 60070 60071 60155 60033 60051 60060 60354 60354 60355 60356 60357 60358 60359 60360 60361 60363 60364 60365 60366 60368 60369 60370 60371 60372 60373 60374 60375 60376 60377 60378 60379 60380 60381 60382 60383 60384 60386 60387 60388 60389 60390 60392 60393 60394 60395 60396 60397 60398 60399 60400 60401 60402 60403 60404 60405 60406 60407 60408 60409 60410 60411 60412 60413 60414 60415 60416 60417 60418 60419 60420 60421 60422 60423 60425 60426 60427 60428 60429 60430 60431 60432 60433 60434 60435 60436 60437 60438 60439 60440 60441 60442 60443 60444 60445 60446 60447 60448 60449 60450 60451 60452 60453 60454 60455 60456 60457 60458 60459 60460 60461 60462 60463 60464 60465 60466 60467 60468 60469 60470 60471 60472 60473 60474 60475 60476 60477 60478 60479 60480 60481 60482 60483 60484 60485 60486 60487 60488 60489 60490 60491 60492 60493 60494 60495 60496 60497 60498 60499 60500 60501 60502 60503 60504 60505 60506 60507 60508 60509 60510 60511 60512 60513 60514 60515 60516 60517 60518 60519 60520 60521 60522 60523 60524 60525 60526 60527 60528 60529 60530 60531 60532 60533 60534 60535 60536 60537 60538 60539 60540 60541 60542 60543 60544 60546 60547 60548 ];
    

%     constants.type.t5 = [60018 60020 60034 60035 60054 60061 60064 60071 60155 60033 60060 60357 60360 60364 60366 60371 60374 60380 60383 60386 60387 60388 60389 60390 60395 60399 60403 60404 60405 60406 60408 60409 60410 60411 60414 60415 60419 60420 60421 60427 60432 60434 60436 60439 60440 60442 60447 60448 60450 60452 60453 60456 60458 60459 60460 60462 60463 60465 60466 60469 60472 60477 60478 60481 60484 60486 60488 60490 60491 60494 60495 60497 60498];
%     constants.type.xbc = [60070 60355 60379 60413 60430 60445 60449 60455 60473 60493];
%     constants.type.t6 = [60017 60031 60036 60356 60358 60359 60363 60365 60369 60370 60372 60375 60381 60382 60392 60394 60398 60400 60401 60407 60416 60423 60425 60426 60428 60431 60441 60443 60444 60451 60464 60467 60468 60471 60475 60476 60479 60483 60487 60496];
%     constants.type.t7 = [60016 60051 60354 60354 60361 60373 60376 60377 60393 60396 60397 60412 60418 60429 60435 60437 60446 60454 60470 60480 60485 60492  ];
%     constants.type.t8 = [60368 60402 60417 60433 60438 60457 60461 60482];
%     constants.type.t9 = [60378 60384 60474];
% 
%     constants.type.t5y = [60020 60054 60064 60071 60155 60060 60357 60364 60366 60371 60386 60387 60388 60389 60399 60403 60405 60406 60411 60415 60421 60427 60436 60442 60450 60456 60459 60462 60465 60469 60477 60478 60484 60490 60494];
%     constants.type.t5z = [60018 60061 60033 60360 60374 60380 60383 60390 60395 60404 60408 60409 60410 60414 60419 60420 60432 60434 60439 60440 60447 60448 60452 60453 60458 60460 60463 60466 60472 60481 60486 60488 60491 60495 60497 60498];

%     constants.type.t5 = [60445 60473 60018 60020 60034 60035 60054 60061 60064 60071 60155 60033   60360 60364 60366 60371 60374 60380 60383 60386 60387 60388 60389 60390 60395 60399 60403 60404 60405 60406 60408  60410 60411 60414 60415  60420 60421  60432 60434 60436 60439 60440 60442 60447 60448 60450 60452 60453  60458 60459 60460 60462  60465 60466 60469   60478 60481 60484  60488 60490 60491 60495 60497 60498];
%     constants.type.xbc = [60070 60355 60379 60413 60430  60449 60455  60493];
%     constants.type.t6 = [60017 60031 60036 60356 60358 60359 60363 60365 60369 60370  60375 60381 60382 60392 60394 60398  60401 60407 60416 60423 60425 60426 60428 60431 60441 60443 60444 60451 60464 60467   60475 60476 60479 60483 60487 60496 60427 60456 60357 60472 60463 60060 60419 60409 60486 60477];
%     constants.type.t7 = [60016 60051 60354 60354 60361 60373 60376 60377 60393 60396 60397 60412 60418 60429 60435 60437 60446 60454 60470 60480 60485 60492  ];
%     constants.type.t8 = [60368 60402 60417 60433 60438 60457 60461 60482];
%     constants.type.rbc = [60378 60384 60474 60468 60400 60471 60372];

    constants.type.xbc = [ 60355 60379 60413 60430  60449 60455 60493 60501 60517  60539 60547 ];
    constants.type.t5w = [60364 60054 60155 60366 60371 60395 60399 60403 60408 60411 60436 60440 60448 60452 60465 60469 60481 60502 60503 60507 60527 60532 60533 60535 60538  ]; %60543 looks non-bipolar?
constants.type.t5l = [60421 60033 60528 60020 60388 60360 60374 60380  60386  60389 60404 60410 60414 60415 60439 60442  60450 60458  60460 60462 60478 60488  60491 60497 60498 60504 60505 60510  60514 60519 60522 60523 60541 60542 ];
constants.type.t5h = [60473  60383 60447 60459 60513 60490 60484  60018 60035 60061 60064 60071 60387 60390 60405 60406 60420  60434 60453 60466 60495  60534 60540 ];
constants.type.t7 = [60016 60051 60354 60354 60361 60373 60376 60377 60393 60396 60397 60412 60418 60429 60435 60437 60446 60454 60470 60480 60485 60492 60508 60509 60525 60529 60531 60546 ];
constants.type.t89 = [60368 60402 60417 60433 60438 60457 60461 60482 60500 ];
constants.type.tRBC = [60017 60365 60369 60372 60378 60384 60392 60400 60463 60468 60471 60474 60476 60506 60518 60521 60544 ];
constants.type.t6 = [60031 60036 60060 60356 60357 60358 60359 60363 60370 60375 60381 60382 60394 60398 60401 60407 60409 60416 60419 60422 60423 60425 60426 60428 60431 60441 60443 60444 60451 60456 60464 60467 60472 60475 60477 60479 60483 60486 60487 60489 60496 60499 60512 60515 60516 60520 60526 60530 60536 60537 60548 ];

    
    
    
    constants.type.coarse_bs = [ 17097 17084 17140 17114];


%     constants.type.g1 = [60034 60070 60427 60494 60511 ];
%     constants.type.g2 = [60355 60379 60413 60430 60449 60455 60493 60501 60517 60539 60547 ];
%     constants.type.g11 = [ 60060 60357 60398 60409 60419 60441 60467 60472 60477 60486 60487 ];
% constants.type.g3 = [ 60054 60155 60366 60371 60395 60399 60403 60408 60411 60436 60440 60445 60448 60452 60465 60469 60481 60502 60503 60507 60524 60527 60532 60533 60535 60538 60543 ];
% constants.type.g4 = [ 60020 60033 60360 60374 60380 60383 60386 60388 60389 60404 60410 60414 60415 60432 60439 60442 60447 60450 60458 60459 60460 60462 60478 60484 60488 60490 60491 60497 60498 60504 60505 60510 60513 60514 60519 60522 60523 60541 60542 ];
% constants.type.g5 = [ 60018 60035 60061 60064 60071 60364 60387 60390 60405 60406 60420 60421 60434 60453 60466 60473 60495 60528 60534 60540 ];
% constants.type.g6 = [ 60016 60051 60354 60354 60361 60373 60376 60377 60393 60396 60397 60412 60418 60429 60435 60437 60446 60454 60470 60480 60485 60492 60508 60509 60525 60529 60531 60546 ];
% constants.type.g7 = [ 60368 60402 60417 60433 60438 60457 60461 60482 60500 ];
% constants.type.g8 = [ 60378 60384 60468 60474 ];
% constants.type.g9 = [ 60017 60031 60036 60356 60358 60359 60363 60370 60372 60375 60382 60392 60394 60400 60401 60407 60416 60423 60428 60444 60464 60479 60483 60489 60496 60499 60506 60516 60518 60526 60530 60537 60544 60548 ];
% constants.type.g10 = [ 60365 60369 60381 60422 60425 60426 60431 60443 60451 60456 60463 60471 60475 60476 60512 60515 60520 60521 60536 ];

% constants.type.g1 = [60034 60070 60427 60494 60511 ];
%     constants.type.g2 = [ 60355 60379 60413 60430 60445 60449 60455 60473 60493 60501 60517 60524 60539 60547 ];
% constants.type.g3 = [ 60054 60155 60366 60371 60395 60399 60403 60408 60411 60436 60440 60448 60452 60465 60469 60481 60502 60503 60507 60527 60532 60533 60535 60538 60543 ];
% constants.type.g4 = [ 60020 60033 60360 60374 60380 60383 60386 60388 60389 60404 60410 60414 60415 60432 60439 60442 60447 60450 60458 60459 60460 60462 60478 60484 60488 60490 60491 60497 60498 60504 60505 60510 60513 60514 60519 60522 60523 60541 60542 ];
% constants.type.g5 = [ 60018 60035 60061 60064 60071 60364 60387 60390 60405 60406 60420 60421 60434 60453 60466 60495 60528 60534 60540 ];
% constants.type.g6 = [ 60016 60051 60354 60354 60361 60373 60376 60377 60393 60396 60397 60412 60418 60429 60435 60437 60446 60454 60470 60480 60485 60492 60508 60509 60525 60529 60531 60546 ];
% constants.type.g7 = [ 60368 60402 60417 60433 60438 60457 60461 60482 60500 ];
% constants.type.g8 = [ 60378 60384 60468 60474 ];
% constants.type.g9 = [ 60017 60031 60036 60356 60358 60359 60363 60370 60372 60375 60382 60392 60394 60400 60401 60407 60416 60423 60428 60444 60464 60479 60483 60489 60496 60499 60506 60516 60518 60526 60530 60537 60544 60548 ];
% constants.type.g10 = [ 60365 60369 60381 60422 60425 60426 60431 60443 60451 60456 60463 60471 60475 60476 60512 60515 60520 60521 60536 ];
% constants.type.g11 = [ 60060 60357 60398 60409 60419 60441 60467 60472 60477 60486 60487 ];





    constants.type.off_sac = [70014 70016 70023 70024 70025 70026 70027 70030 70031 70032 70033 70034 70035 70048 70050 ...
        70066 70068 70076 70077 70079 70080 70081 70082 70083 70084 70085 70086 70087 70088 70089 70090 70093 70094 70095 ...
        70096 70097 70098 70099 70100 70101 70102 70104 70105 70106 70108 70109 70110 70111 70112 70113 70114 70115 70116 ...
        70117 70118 70119 70120 70121 70122 70123 70124 70125 70126 70127 70128 70129 70130 70131 70132 70133 70134 70135 ...
        70136 70137 70138 70139 70140 70141 70142 70143 70144 70145 70146 70147 70148 70149 70150 70151 70152 70154 70155 ...
        70156 70158];


    


    constants.type.sure_off_sac = [70014       70016       70023       70024       70025 ...
       70026       70027       70030       70031       70032 ...
       70033       70034       70035       70048       70050 ...
       70066       70068       70076       70077       70079 ...
       70080       70081       70082       70083       70084 ...
       70085       70086       70087       70088       70089 ...
       70090       70093       70095       70096       70099 ...
       70100       70102       70106       70108       70109 ...
       70110       70111       70112       70113       70114 ...
       70115       70116       70117       70118       70119 ...
       70120       70121       70122       70123       70124 ...
       70125       70126       70127       70128       70129 ...
       70130       70131       70133       70134       70137 ...
       70138       70141       70145       70146       70147 ...
       70148       70149       70150       70151       70152 ...
70154       70155       70156       70158];

    constants.type.sure_t3a = [60024 60028 60048 60049 60050 60066 60068 60072 60075 60076 60085 60092 60093 60096 60108 60123 60134 60136 60146 60166 60181];
    constants.type.sure_t4 = [60059 60083 60088 60089 60116 60163 60186];


    
    constants.type.C1 = [           60008	       60019	       60026	       60027	       60032	       60052	       60055	       60078	       60079	       60099	       60105	       60109	       60110	       60111	       60114	       60118	       60129	       60132	       60139	       60147	       60150	       60158	       60162	       60164	       60170	       60177	       60184	       60187	       60189	       60196	       60203	       60216	       60218	       60230	       60248	       60272	       60277	       60279	       60187	       60217];
    constants.type.C2 = [60142	       60161	       60194	       60195	       60204	       60227	       60012	       60039	       60104	       60157	       60169	       60182	       60214	       60221	       60226	       60236	       60257	       60269];
    constants.type.C3 = [60188	       60212	       60213	       60001	       60002	       60003	       60004	       60006	       60009	       60010	       60011	       60013	       60021	       60022	       60023	       60025	       60029	       60037	       60038	       60040	       60041	       60042	       60043	       60044	       60046	       60080	       60097	       60101	       60102	       60103	       60112	       60120	       60124	       60130	       60133	       60135	       60138	       60140	       60141	       60149	       60159	       60167	       60198	       60208	       60209	       60210	       60219	       60223	       60224	       60235	       60241	       60242	       60245	       60254	       60261	       60263	       60266	       60267	       60268	       60270	       60274	       60285	       60286];

    
% constants.type.on_sac = [70000+[29] 1070000+[28 161 162 163 164 168:172 174 176 178:189 191:209]];
% constants.type.on_sac = [70029 70028 70161 70162 70163 70164 70168 70169 70170 70171 70172 70174 70176 70178 70179 70180 70181 70182 70183 70184 70185 70186 70187 70188 70189 70191 70192 70193 70194 70195 70196 70197 70198 70199 70200 70201 70202 70203 70204 70205 70206 70207 70208 70209 70211 70212 70213 70214 70215 70216 70217 70218 70219 70220 70221 70222 70223 70224 70225 70227 70228];% 70232 70233 70234 70235 70236 70237 70238 70239 70240 70241 70242 70243 70244];

constants.type.on_sac = [70028 70029 70161 70162 70163 70164 70168 70169 70170 70171 70172 70174 70176 70179 70180 70181 70182 70183 70184 70185 70186 70187 70188 70189 70191 70192 70193 70194 70195 70196 70197 70198 70199 70200 70201 70202 70203 70204 70205 70206 70207 70208 70209 70211 70212 70213 70214 70215 70216 70217 70218 70219 70220 70221 70222 70223 70224 70225 70227 70228 70229 70230 70231 70232 70233 70234 70235 70236 70237 70238 70239 70240 70241 70242 70243 70244 ]; %70178 omitted for merger

% 
%      
%      [70014       70016       70048       70050       70066       70068 ...
%          70076       70077       70079       70080       70081       70082 ...
%          70083       70084       70085       70086       70087       70088 ...
%          70089       70090       70093       70094       70095       70096 ...
%          70097       70098       70099       70100       70101       70102 ...
%         70104       70105       70106       70108       70023       70024 ...
%         70025       70026       70027       70030       70031       70032 ...
%          70033       70034       70035];
     
     constants.type.forward_sac = [70023 70024 70025 70026 70030 70031 ...
         70032 70033 70034 70035 70048 70050 70066 70076 70077 70079 ...
         70080 70084 70087 70090 70095 70106];
     constants.type.backward_sac = [70023, 70024, 70025, 70026, 70031, ...
         70032, 70033, 70035, 70048, 70050, 70066, 70068, 70076, 70077, ...
         70079, 70084, 70087, 70095, 70106, 70109, 70110, 70111, 70112, ...
         70113, 70114, 70116, 70117, 70118, 70119, 70121, 70122, 70123, ...
         70124, 70125, 70127, 70128, 70130];
     constants.type.minij = [10010, 10017, 15066, 15018, 17027, 17105 17177];
     constants.type.j = [17028 17060 17075];
    
    
end