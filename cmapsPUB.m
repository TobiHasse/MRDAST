function my_colors = cmapsPUB(my_map,flp)
% Purpose:  This function will return colormaps used in various figures in
%           Tobias Hasse's dissertation 

% Inputs:   my_map          the name of my colormap as a string
%           flp             boolean whether to flip (reverse) colormap
% Actions:                  set the colormap on gcf
% Outputs:  my_colors       (optional) rgb array of colormap data
% Author:   Tobias Hasse    tobiack@udel.edu
% Date:     June 2021

% colormaps
%This file is to preserve the colormaps for different figures
% I think a switch & case could be used to return the appropriate colormap
%possibly each map could be listed in its own function that could be called
%separately or from this main function

switch my_map
    case 'age'
%         colormap(cmap_age)
        my_colors = cmap_age;
    case 'topo'
%         colormap(cmap_topo)
        my_colors = cmap_topo;
    case 'xsec'
        my_colors = cmap_xsec;
    case 'median_stripe'
        my_colors = cmap_topo;
        my_colors = cmap_median(my_colors);
    otherwise
        fprintf(strcat('colormap choice failed, try:',...
            ' age, topo, xsec, or median_stripe\n'))
end
if flp
    my_colors = flipud(my_colors);
end
colormap(my_colors)

end


function cmap = cmap_median(cmap) % replace the middle with cyan stripe
% make stripe width depend on colormap width
id_mid = round(length(cmap)/2);
idx =ceil(length(cmap)/100);
% id_mid-idx;
% id_mid+idx;
% keyboard
cmap(id_mid-idx+1 : id_mid+idx,:) = repmat([ 0 1 1 ], 2*idx,1);
% cdata2(126:131,:)=[0 1 1;0 1 1;0 1 1;0 1 1;0 1 1;0 1 1]; % for 256 colors

end


function cmap = cmap_age()
% age colormap
cmap= [  0         0    1.0000;
    0.8929    0.6945    0.3969;
    0.8787    0.6835    0.3906;
    0.8646    0.6724    0.3843;
    0.8504    0.6614    0.3780;
    0.8362    0.6504    0.3717;
    0.8220    0.6394    0.3654;
    0.8079    0.6283    0.3591;
    0.7937    0.6173    0.3528;
    0.7654    0.5953    0.3402;
    0.7512    0.5843    0.3339;
    0.7370    0.5732    0.3276;
    0.7228    0.5622    0.3213;
    0.7087    0.5512    0.3150;
    0.6945    0.5402    0.3087;
    0.6803    0.5291    0.3024;
    0.6661    0.5181    0.2961;
    0.6520    0.5071    0.2898;
    0.6236    0.4850    0.2772;
    0.6094    0.4740    0.2709;
    0.5953    0.4630    0.2646;
    0.5811    0.4520    0.2583;
    0.5669    0.4409    0.2520;
    0.5528    0.4299    0.2457;
    0.5386    0.4189    0.2394;
    0.5244    0.4079    0.2331;
    0.5102    0.3969    0.2268;
    0.4819    0.3748    0.2142;
    0.4677    0.3638    0.2079;
    0.4535    0.3528    0.2016;
    0.4394    0.3417    0.1953;
    0.4252    0.3307    0.1890;
    0.4110    0.3197    0.1827;
    0.3969    0.3087    0.1764;
    0.3827    0.2976    0.1701;
    0.3685    0.2866    0.1638;
    0.3402    0.2646    0.1512;
    0.3260    0.2535    0.1449;
    0.3118    0.2425    0.1386;
    0.2976    0.2315    0.1323;
    0.2835    0.2205    0.1260;
    0.2693    0.2094    0.1197;
    0.2551    0.1984    0.1134;
    0.2409    0.1874    0.1071;
    0.2268    0.1764    0.1008;
    0.2126    0.1654    0.0945;
    0.1984    0.1543    0.0882;
    0.1843    0.1433    0.0819;
    0.1701    0.1323    0.0756;
    0.1559    0.1213    0.0693;
    0.1417    0.1102    0.0630;
    0.1276    0.0992    0.0567;
    0.1134    0.0882    0.0504;
    0.0992    0.0772    0.0441;
    0.0850    0.0661    0.0378;
    0.0709    0.0551    0.0315;
    0.0567    0.0441    0.0252;
    0.0425    0.0331    0.0189;
    0.0283    0.0220    0.0126;
    0.0142    0.0110    0.0063;
         0         0         0;
    0.0069    0.0069    0.0121;
    0.0207    0.0207    0.0311;
    0.0344    0.0344    0.0501;
    0.0482    0.0482    0.0691;
    0.0620    0.0620    0.0880;
    0.0758    0.0758    0.1070;
    0.0896    0.0896    0.1260;
    0.1033    0.1033    0.1450;
    0.1171    0.1171    0.1640;
    0.1309    0.1309    0.1830;
    0.1447    0.1447    0.2020;
    0.1585    0.1585    0.2210;
    0.1722    0.1722    0.2400;
    0.1860    0.1860    0.2589;
    0.1998    0.1998    0.2779;
    0.2136    0.2136    0.2969;
    0.2274    0.2274    0.3159;
    0.2411    0.2411    0.3349;
    0.2549    0.2549    0.3539;
    0.2687    0.2687    0.3729;
    0.2825    0.2825    0.3919;
    0.2963    0.2963    0.4108;
    0.3100    0.3100    0.4298;
    0.3238    0.3238    0.4488;
    0.3376    0.3428    0.4626;
    0.3514    0.3618    0.4764;
    0.3652    0.3808    0.4902;
    0.3789    0.3998    0.5039;
    0.3927    0.4188    0.5177;
    0.4065    0.4377    0.5315;
    0.4203    0.4567    0.5453;
    0.4341    0.4757    0.5591;
    0.4478    0.4947    0.5728;
    0.4616    0.5137    0.5866;
    0.4754    0.5327    0.6004;
    0.4892    0.5517    0.6142;
    0.5030    0.5707    0.6280;
    0.5167    0.5896    0.6417;
    0.5305    0.6086    0.6555;
    0.5443    0.6276    0.6693;
    0.5581    0.6466    0.6831;
    0.5719    0.6656    0.6969;
    0.5856    0.6846    0.7106;
    0.5994    0.7036    0.7244;
    0.6132    0.7226    0.7382;
    0.6270    0.7416    0.7520;
    0.6407    0.7605    0.7657;
    0.6545    0.7795    0.7795;
    0.6761    0.7933    0.7933;
    0.6977    0.8071    0.8071;
    0.7193    0.8209    0.8209;
    [0.7409    0.8346    0.8346];
    [0.7625    0.8484    0.8484];
    [0.7841    0.8622    0.8622];
    [0.8057    0.8760    0.8760];
    [0.8273    0.8898    0.8898];
    [0.8489    0.9035    0.9035];
    [0.8704    0.9173    0.9173];
    [0.8900    0.9300    0.9300];
    [0.8922    0.9256    0.9456];
    [0.8944    0.9211    0.9611];
    [0.8967    0.9167    0.9767];
    [0.8989    0.9122    0.9922];
    [0.8959    0.8976    1.0000];
    [0.8878    0.8727    1.0000];
    [0.8796    0.8478    1.0000];
    [0.8714    0.8229    1.0000];
    [0.8633    0.7980    1.0000];
    0.8551    0.7731    1.0000;
    0.8469    0.7482    1.0000;
    0.8388    0.7233    1.0000;
    0.8306    0.6984    1.0000;
    0.8224    0.6735    1.0000;
    0.8143    0.6486    1.0000;
    0.8061    0.6237    1.0000;
    0.7980    0.5988    1.0000;
    0.7898    0.5739    1.0000;
    0.7816    0.5490    1.0000;
    0.7735    0.5241    1.0000;
    0.7653    0.4992    1.0000;
    0.7571    0.4743    1.0000;
    0.7490    0.4494    1.0000;
    0.7408    0.4245    1.0000;
    0.7327    0.3996    1.0000;
    0.7245    0.3747    1.0000;
    0.7163    0.3498    1.0000;
    0.7082    0.3249    1.0000;
    0.7000    0.3000    1.0000;
    0.6667    0.3333    1.0000;
    0.6349    0.3651    1.0000;
    0.6032    0.3968    1.0000;
    0.5714    0.4286    1.0000;
    0.5397    0.4603    1.0000;
    0.5079    0.4921    1.0000;
    0.4762    0.5238    1.0000;
    0.4444    0.5556    1.0000;
    0.4127    0.5873    1.0000;
    0.3810    0.6190    1.0000;
    0.3492    0.6508    1.0000;
    0.3175    0.6825    1.0000;
    0.2857    0.7143    1.0000;
    0.2540    0.7460    1.0000;
    0.2222    0.7778    1.0000;
    0.1905    0.8095    1.0000;
    0.1587    0.8413    1.0000;
    0.1270    0.8730    1.0000;
    0.0952    0.9048    1.0000;
    0.0635    0.9365    1.0000;
    0.0317    0.9683    1.0000;
         0    1.0000    1.0000;
         0    0.9792    0.9750;
         0    0.9375    0.9250;
         0    0.8958    0.8750;
         0    0.8542    0.8250;
         0    0.8125    0.7750;
         0    0.7708    0.7250;
         0    0.7292    0.6750;
         0    0.6875    0.6250;
         0    0.6458    0.5750;
         0    0.6042    0.5250;
         0    0.5625    0.4750;
         0    0.5208    0.4250;
    0.0159    0.5079    0.4000;
    0.0476    0.5238    0.4000;
    0.0794    0.5397    0.4000;
    0.1111    0.5556    0.4000;
    0.1429    0.5714    0.4000;
    0.1746    0.5873    0.4000;
    0.2063    0.6032    0.4000;
    0.2381    0.6190    0.4000;
    0.2698    0.6349    0.4000;
    0.3016    0.6508    0.4000;
    0.3333    0.6667    0.4000;
    0.3651    0.6825    0.4000;
    0.3968    0.6984    0.4000;
    0.4286    0.7143    0.4000;
    0.4603    0.7302    0.4000;
    0.4921    0.7460    0.4000;
    0.5238    0.7619    0.4000;
    0.5556    0.7778    0.4000;
    0.5873    0.7937    0.4000;
    0.6190    0.8095    0.4000;
    0.6508    0.8254    0.4000;
    0.6825    0.8413    0.4000;
    0.7143    0.8571    0.4000;
    0.7460    0.8730    0.4000;
    0.7778    0.8889    0.4000;
    0.8095    0.9048    0.4000;
    0.8413    0.9206    0.4000;
    0.8730    0.9365    0.4000;
    0.9048    0.9524    0.4000;
    0.9365    0.9683    0.4000;
    0.9683    0.9841    0.4000;
    1.0000    1.0000    0.4000;
    1.0000    0.9920    0.0003;
    1.0000    0.9759    0.0023;
    1.0000    0.9595    0.0063;
    1.0000    0.9428    0.0123;
    1.0000    0.9258    0.0204;
    1.0000    0.9085    0.0305;
    1.0000    0.8909    0.0426;
    1.0000    0.8729    0.0567;
    1.0000    0.8545    0.0728;
    1.0000    0.8357    0.0910;
    1.0000    0.8165    0.1111;
    1.0000    0.7968    0.1333;
    1.0000    0.7766    0.1575;
    1.0000    0.7559    0.1837;
    1.0000    0.7346    0.2119;
    1.0000    0.7127    0.2421;
    1.0000    0.6901    0.2744;
    1.0000    0.6667    0.3086;
    1.0000    0.6424    0.3449;
    1.0000    0.6172    0.3832;
    1.0000    0.5909    0.4235;
    1.0000    0.5634    0.4659;
    1.0000    0.5345    0.5102;
    1.0000    0.5040    0.5566;
    1.0000    0.4714    0.6049;
    1.0000    0.4364    0.6553;
    1.0000    0.3984    0.7077;
    1.0000    0.3563    0.7622;
    1.0000    0.3086    0.8186;
    1.0000    0.2520    0.8770;
    1.0000    0.1782    0.9375;
    1.0000         0    1.0000;
    1.0000         0    0.5000;
    1.0000    0.1969    0.1969;
    1.0000    0.3807    0.3807;
    1.0000    0.5173    0.5173;
    1.0000    0.6330    0.6330;
    1.0000    0.7360    0.7360;
    1.0000    0.8302    0.8302;
    1.0000    0.9177    0.9177;
    1.0000    1.0000    1.0000];

% make transition from light gray to lavendar a bit darker
my_scalar=ones(length(cmap),1);
grayer = 2:2:16;
grayer = [grayer,18,fliplr(grayer)];
grayer(grayer>14)=14;
my_scalar(113:129) = my_scalar(113:129) + grayer'/100;

% cmap_grayer = bsxfun(@rdivide,cmap,my_scalar);
cmap = bsxfun(@rdivide,cmap,my_scalar);
% colormap(cmap_grayer)
end


function cmap = cmap_topo()
cmap = [ 0         0    1.0000;
    0.8100    0.7000    0.6200;
    0.8059    0.6933    0.6116;
    0.8019    0.6866    0.6031;
    0.7978    0.6798    0.5947;
    0.7938    0.6731    0.5863;
    0.7897    0.6664    0.5778;
    0.7856    0.6597    0.5694;
    0.7816    0.6530    0.5609;
    0.7775    0.6462    0.5525;
    0.7734    0.6395    0.5441;
    0.7694    0.6328    0.5356;
    0.7653    0.6261    0.5272;
    0.7613    0.6194    0.5188;
    0.7572    0.6127    0.5103;
    0.7531    0.6059    0.5019;
    0.7491    0.5992    0.4934;
    0.7450    0.5925    0.4850;
    0.7409    0.5858    0.4766;
    0.7369    0.5791    0.4681;
    0.7328    0.5723    0.4597;
    0.7288    0.5656    0.4513;
    0.7247    0.5589    0.4428;
    0.7206    0.5522    0.4344;
    0.7166    0.5455    0.4259;
    0.7125    0.5388    0.4175;
    0.7084    0.5320    0.4091;
    0.7044    0.5253    0.4006;
    0.7003    0.5186    0.3922;
    0.6963    0.5119    0.3838;
    0.6922    0.5052    0.3753;
    0.6881    0.4984    0.3669;
    0.6841    0.4917    0.3584;
    0.6800    0.4850    0.3500;
    0.6759    0.4783    0.3416;
    0.6719    0.4716    0.3331;
    0.6678    0.4648    0.3247;
    0.6638    0.4581    0.3162;
    0.6597    0.4514    0.3078;
    0.6556    0.4447    0.2994;
    0.6516    0.4380    0.2909;
    0.6475    0.4312    0.2825;
    0.6434    0.4245    0.2741;
    0.6394    0.4178    0.2656;
    0.6353    0.4111    0.2572;
    0.6313    0.4044    0.2488;
    0.6272    0.3977    0.2403;
    0.6231    0.3909    0.2319;
    0.6191    0.3842    0.2234;
    0.6150    0.3775    0.2150;
    0.6109    0.3708    0.2066;
    0.6069    0.3641    0.1981;
    0.6028    0.3573    0.1897;
    0.5988    0.3506    0.1812;
    0.5947    0.3439    0.1728;
    0.5906    0.3372    0.1644;
    0.5866    0.3305    0.1559;
    0.5825    0.3237    0.1475;
    0.5784    0.3170    0.1391;
    0.5744    0.3103    0.1306;
    0.5703    0.3036    0.1222;
    0.5663    0.2969    0.1137;
    0.5622    0.2902    0.1053;
    0.5581    0.2834    0.0969;
    0.5541    0.2767    0.0884;
    0.5500    0.2700    0.0800;
    0.5570    0.2658    0.0788;
    0.5641    0.2616    0.0775;
    0.5711    0.2573    0.0762;
    0.5781    0.2531    0.0750;
    0.5852    0.2489    0.0737;
    0.5922    0.2447    0.0725;
    0.5992    0.2405    0.0713;
    0.6063    0.2363    0.0700;
    0.6133    0.2320    0.0688;
    0.6203    0.2278    0.0675;
    0.6273    0.2236    0.0663;
    0.6344    0.2194    0.0650;
    0.6414    0.2152    0.0638;
    0.6484    0.2109    0.0625;
    0.6555    0.2067    0.0612;
    0.6625    0.2025    0.0600;
    0.6695    0.1983    0.0588;
    0.6766    0.1941    0.0575;
    0.6836    0.1898    0.0563;
    0.6906    0.1856    0.0550;
    0.6977    0.1814    0.0537;
    0.7047    0.1772    0.0525;
    0.7117    0.1730    0.0513;
    0.7188    0.1688    0.0500;
    0.7258    0.1645    0.0488;
    0.7328    0.1603    0.0475;
    0.7398    0.1561    0.0462;
    0.7469    0.1519    0.0450;
    0.7539    0.1477    0.0438;
    0.7609    0.1434    0.0425;
    0.7680    0.1392    0.0413;
    0.7750    0.1350    0.0400;
    0.7820    0.1308    0.0388;
    0.7891    0.1266    0.0375;
    0.7961    0.1223    0.0362;
    0.8031    0.1181    0.0350;
    0.8102    0.1139    0.0338;
    0.8172    0.1097    0.0325;
    0.8242    0.1055    0.0313;
    0.8313    0.1013    0.0300;
    0.8383    0.0970    0.0288;
    0.8453    0.0928    0.0275;
    0.8523    0.0886    0.0262;
    0.8594    0.0844    0.0250;
    0.8664    0.0802    0.0238;
    0.8734    0.0759    0.0225;
    0.8805    0.0717    0.0213;
    0.8875    0.0675    0.0200;
    0.8945    0.0633    0.0187;
    0.9016    0.0591    0.0175;
    0.9086    0.0548    0.0163;
    0.9156    0.0506    0.0150;
    0.9227    0.0464    0.0138;
    0.9297    0.0422    0.0125;
    0.9367    0.0380    0.0113;
    0.9437    0.0338    0.0100;
    0.9508    0.0295    0.0088;
    0.9578    0.0253    0.0075;
    0.9648    0.0211    0.0063;
    0.9719    0.0169    0.0050;
    0.9789    0.0127    0.0037;
    0.9859    0.0084    0.0025;
    0.9930    0.0042    0.0013;
    1.0000         0         0;
    1.0000    0.0149         0;
    1.0000    0.0299         0;
    1.0000    0.0448         0;
    1.0000    0.0597         0;
    1.0000    0.0746         0;
    1.0000    0.0896         0;
    1.0000    0.1045         0;
    1.0000    0.1194         0;
    1.0000    0.1343         0;
    1.0000    0.1493         0;
    1.0000    0.1642         0;
    1.0000    0.1791         0;
    1.0000    0.1940         0;
    1.0000    0.2090         0;
    1.0000    0.2239         0;
    1.0000    0.2388         0;
    1.0000    0.2537         0;
    1.0000    0.2687         0;
    1.0000    0.2836         0;
    1.0000    0.2985         0;
    1.0000    0.3134         0;
    1.0000    0.3284         0;
    1.0000    0.3433         0;
    1.0000    0.3582         0;
    1.0000    0.3731         0;
    1.0000    0.3881         0;
    1.0000    0.4030         0;
    1.0000    0.4179         0;
    1.0000    0.4328         0;
    1.0000    0.4478         0;
    1.0000    0.4627         0;
    1.0000    0.4776         0;
    1.0000    0.4925         0;
    1.0000    0.5075         0;
    1.0000    0.5224         0;
    1.0000    0.5373         0;
    1.0000    0.5522         0;
    1.0000    0.5672         0;
    1.0000    0.5821         0;
    1.0000    0.5970         0;
    1.0000    0.6119         0;
    1.0000    0.6269         0;
    1.0000    0.6418         0;
    1.0000    0.6567         0;
    1.0000    0.6716         0;
    1.0000    0.6866         0;
    1.0000    0.7015         0;
    1.0000    0.7164         0;
    1.0000    0.7313         0;
    1.0000    0.7463         0;
    1.0000    0.7612         0;
    1.0000    0.7761         0;
    1.0000    0.7910         0;
    1.0000    0.8060         0;
    1.0000    0.8209         0;
    1.0000    0.8358         0;
    1.0000    0.8507         0;
    1.0000    0.8657         0;
    1.0000    0.8806         0;
    1.0000    0.8955         0;
    1.0000    0.9104         0;
    1.0000    0.9254         0;
    1.0000    0.9403         0;
    1.0000    0.9552         0;
    1.0000    0.9701         0;
    1.0000    0.9851         0;
    1.0000    1.0000         0;
    0.9167    1.0000         0;
    0.8333    1.0000         0;
    0.7500    1.0000         0;
    0.6667    1.0000         0;
    0.5833    1.0000         0;
    0.5000    1.0000         0;
    0.4167    1.0000         0;
    0.3333    1.0000         0;
    0.2500    1.0000         0;
    0.1667    1.0000         0;
    0.0833    1.0000         0;
         0    1.0000         0;
         0    1.0000    0.0909;
         0    1.0000    0.1818;
         0    1.0000    0.2727;
         0    1.0000    0.3636;
         0    1.0000    0.4545;
         0    1.0000    0.5455;
         0    1.0000    0.6364;
         0    1.0000    0.7273;
         0    1.0000    0.8182;
         0    1.0000    0.9091;
         0    1.0000    1.0000;
         0    0.8889    1.0000;
         0    0.7778    1.0000;
         0    0.6667    1.0000;
         0    0.5556    1.0000;
         0    0.4444    1.0000;
         0    0.3333    1.0000;
         0    0.2222    1.0000;
         0    0.1111    1.0000;
         0         0    1.0000;
    0.0909         0    1.0000;
    0.1818         0    1.0000;
    0.2727         0    1.0000;
    0.3636         0    1.0000;
    0.4545         0    1.0000;
    0.5455         0    1.0000;
    0.6364         0    1.0000;
    0.7273         0    1.0000;
    0.8182         0    1.0000;
    0.8600         0    1.0000;  %    0.9091         0    1.0000;
    1.0000       .12    1.0000;  %    1.0000         0    1.0000;
    0.9400    0.2100    0.9400;
    0.8525    0.2138    0.8525;
    0.7650    0.2175    0.7650;
    0.6775    0.2213    0.6775;
    0.5900    0.2250    0.5900;
    0.5025    0.2288    0.5025;
    0.4150    0.2325    0.4150;
    0.3275    0.2362    0.3275;
    0.2400    0.2400    0.2400;
    0.3486    0.3486    0.3486;
    0.4571    0.4571    0.4571;
    0.5657    0.5657    0.5657;
    0.6743    0.6743    0.6743;
    0.7829    0.7829    0.7829;
    0.8914    0.8914    0.8914;
    1.0000    1.0000    1.0000];
    

end


function cmap = cmap_xsec()
cmap = [1	0.900653595	0.933986928;
1       	0.80130719	0.867973856;
1           0.701960784	0.801960784;
1        	0.602614379	0.735947712;
1           0.503267974	0.669934641;
1       	0.403921569	0.603921569;
0.92    	0.363921569	0.602352941;
0.84        0.323921569	0.600784314;
0.76    	0.283921569	0.599215686;
0.68    	0.243921569	0.597647059;
0.6         0.203921569	0.596078431;
0.514285714	0.317647059	0.653781513;
0.428571429	0.431372549	0.711484594;
0.342857143	0.545098039	0.769187675;
0.257142857	0.658823529	0.826890756;
0.171428571	0.77254902	0.884593838;
0.085714286	0.88627451	0.942296919;
0       	1       	1;
0.111111111	0.954684096	0.912854031;
0.222222222	0.909368192	0.825708061;
0.333333333	0.864052288	0.738562092;
0.444444444	0.818736383	0.651416122;
0.555555556	0.773420479	0.564270153;
0.666666667	0.728104575	0.477124183;
0.777777778	0.682788671	0.389978214;
0.888888889	0.637472767	0.302832244;
1       	0.592156863	0.215686275;
1       	0.637472767	0.191721133;
1       	0.682788671	0.167755991;
1       	0.728104575	0.14379085;
1       	0.773420479	0.119825708;
1       	0.818736383	0.095860566;
1       	0.864052288	0.071895425;
1       	0.909368192	0.047930283;
1       	0.954684096	0.023965142;
1       	1       	0;
0.97254902	0.928104575	0.009150327;
0.945098039	0.85620915	0.018300654;
0.917647059	0.784313725	0.02745098;
0.837472767	0.720261438	0.024400871;
0.757298475	0.65620915	0.021350763;
0.677124183	0.592156863	0.018300654;
0.596949891	0.528104575	0.015250545;
0.516775599	0.464052288	0.012200436;
0.436601307	0.4     	0.009150327;
0.356427015	0.335947712	0.006100218;
0.276252723	0.271895425	0.003050109;
0.196078431	0.207843137	0;
0.098039216	0.353921569	0.2;
0       	0.5     	0.4;
0.089215686	0.543627451	0.391176471;
0.178431373	0.587254902	0.382352941;
0.267647059	0.630882353	0.373529412;
0.356862745	0.674509804	0.364705882;
0.047058824	0.701960784	0.015686275;
0.02745098	0.698039216	0.490196078;
0.11372549	0.654901961	0.543137255;
0.2     	0.611764706	0.596078431;
0.177777778	0.54379085	0.640958606;
0.155555556	0.475816993	0.68583878;
0.133333333	0.407843137	0.730718954;
0.111111111	0.339869281	0.775599129;
0.088888889	0.271895425	0.820479303;
0.066666667	0.203921569	0.865359477;
0.044444444	0.135947712	0.910239651;
0.022222222	0.067973856	0.955119826;
0       	0       	1;
0       	0.006535948	0.797385621;
0       	0.013071895	0.594771242;
0       	0.019607843	0.392156863;
0.086834734	0.101960784	0.478991597;
0.173669468	0.184313725	0.565826331;
0.260504202	0.266666667	0.652661064;
0.347338936	0.349019608	0.739495798;
0.434173669	0.431372549	0.826330532;
0.521008403	0.51372549	0.913165266;
0.607843137	0.596078431	1;
0.606722689	0.511484594	0.88627451;
0.605602241	0.426890756	0.77254902;
0.604481793	0.342296919	0.658823529;
0.603361345	0.257703081	0.545098039;
0.602240896	0.173109244	0.431372549;
0.601120448	0.088515406	0.317647059;
0.6     	0.003921569	0.203921569;
0.733333333	0.002614379	0.135947712;
0.866666667	0.00130719	0.067973856;
1       	0       	0;
1       	0.201960784	0.301960784;
1       	0.302941176	0.452941176;
1       	0.403921569	0.603921569;
0.94379085	0.360348584	0.550762527;
0.887581699	0.316775599	0.497603486;
0.831372549	0.273202614	0.444444444;
0.775163399	0.22962963	0.391285403;
0.718954248	0.186056645	0.338126362;
0.662745098	0.14248366	0.28496732;
0.606535948	0.098910675	0.231808279;
0.550326797	0.055337691	0.178649237;
0.494117647	0.011764706	0.125490196;
0.459183007	0.073300654	0.189281046;
0.424248366	0.134836601	0.253071895;
0.389313725	0.196372549	0.316862745;
0.354379085	0.257908497	0.380653595;
0.319444444	0.319444444	0.444444444;
0.375   	0.395833333	0.5;
0.430555556	0.472222222	0.555555556;
0.486111111	0.548611111	0.611111111;
0.541666667	0.625   	0.666666667;
0.597222222	0.701388889	0.722222222;
0.652777778	0.777777778	0.777777778;
0.739583333	0.833333333	0.833333333;
0.826388889	0.888888889	0.888888889;
0.913194444	0.944444444	0.944444444;
1	        1	        1 ];

end