/*
 * loexp - find expectile curves of local polynomial regression 
 *
 *  Copyright (C) 1999-2008, Pratyaksha J. Wirapati
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 51 Franklin Street
 *  Fifth Floor, Boston, MA 02110-1301  USA.

*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "loexp.h"

/*
 * Utilities
 */ 

// allocate double vector
// 
static inline double*
alloc(int n)
{
  return (double*)malloc(sizeof(double)*n);
}

// set the value of a double vector to the constant c
//
static void
set_d (int n, double c, double *x)
{
  for(int i = 0; i < n; i++ )
    x[i] = c;
}

// Recursive Gaussian filter
//
// ref: Young & van Vliet (1995) "Recursive implementation of the
// Gaussian filter" Signal Processing 44:139-151
//
// Note: the filter coeffs according to the paper above does not
// work in the context regression weights because of small negative
// values in the impulse response (causing
// the normal eq. for local polynomial to have non-positive definite
// X'WX). We use tabulated values of coefficients (see below)
//

// 
// filter the columns of a densely packed matrix
// 
static void
recgauss_filter_col (int n, int m, double *y, int K, double *a)
{
  double *yi = y + K*m;
  for(int i = K; i < n; i++, yi += m )
    {
    for(int j = 0; j < m; j++ )
      yi[j] *= a[0];
    double *yk = yi - m;
    for(int k = 1; k <= K; k++, yk -= m )
      for(int j = 0; j < m; j++ )
        yi[j] -= a[k]*yk[j];
    }
  yi = y + (n-K-1)*m;
  for(int i = n-K-1; i >= K; i--, yi -= m )
    {
    for(int j = 0; j < m; j++ )
      yi[j] *= a[0];
    double *yk = yi + m;
    for(int k = 1; k <= K; k++, yk += m )
      for(int j = 0; j < m; j++ )
        yi[j] -= a[k]*yk[j];
    } 
}

// To ensure that the Gaussian impulse response is non-negative,
// the filter coefficients were found using generic optimization,
// discarding solutions with any negative response.
// 
// In practice, discrete sigma is enough and it's easier to
// use tabulated filter coefficients. Only kernel width up to
// RECGAUSS_MAX_SIGMA is supported. It is sufficient for most applications,
// unless the curve is meant to insensitive to very wide bumps.
// 
// This limit has no significance. It's chosen because it gets more
// difficult to find the coefficients for larger sigma.
//
//
#define RECGAUSS_MAX_SIGMA 283
static double
recgauss_3coeff_nonneg_table[RECGAUSS_MAX_SIGMA+1][4] =
{{1, 0, 0, 0, }
,{4.530743479429e-01, -0.830067524045,  0.376336833667, -0.093194961679}
,{2.249561005778e-01, -1.281213160862,  0.618875754286, -0.112706492845}
,{1.027322175680e-01, -1.708603653126,  1.031177211289, -0.219841340595}
,{5.330075904580e-02, -1.979977701297,  1.350847329385, -0.317568869042}
,{3.123259224558e-02, -2.154139806956,  1.579862916112, -0.394490516911}
,{1.998247542936e-02, -2.275114113434,  1.750885285902, -0.455788697039}
,{1.328137196501e-02, -2.371301168399,  1.894645464894, -0.510062924530}
,{9.218755377370e-03, -2.446853682949,  2.012314130361, -0.556241692035}
,{6.754623239290e-03, -2.502242766814,  2.100702634388, -0.591705244335}
,{5.367114282492e-03, -2.532217382062,  2.147731629047, -0.610147132702}
,{4.156606944747e-03, -2.570222573266,  2.210674652405, -0.636295472194}
,{3.254181896531e-03, -2.605114511419,  2.269743233276, -0.661374539961}
,{2.470606611620e-03, -2.646774952229,  2.342430155623, -0.693184596782}
,{1.978299693897e-03, -2.673304198183,  2.388635972328, -0.713353474452}
,{1.626041545312e-03, -2.694437556777,  2.425707149046, -0.729643550723}
,{1.354284932138e-03, -2.712881691522,  2.458333207055, -0.744097230600}
,{1.160791332206e-03, -2.726672614749,  2.482770037911, -0.754936631829}
,{1.067496270882e-03, -2.728566373788,  2.485338976205, -0.755705106146}
,{9.060093989016e-04, -2.743444369488,  2.512139416031, -0.767789037144}
,{7.160291071272e-04, -2.768600097549,  2.558436110664, -0.789119984007}
,{6.256148116213e-04, -2.778713499832,  2.576836803922, -0.797497689279}
,{5.424653049162e-04, -2.789404242513,  2.596433030093, -0.806486322275}
,{4.769960286316e-04, -2.798353279132,  2.612875336743, -0.814045061582}
,{4.232496515650e-04, -2.806243273796,  2.627414656333, -0.820748132885}
,{3.751215136867e-04, -2.814026882654,  2.641819489052, -0.827417484884}
,{3.353441345497e-04, -2.820888237217,  2.654548040171, -0.833324458819}
,{3.069396653252e-04, -2.825391423925,  2.662853760140, -0.837155396550}
,{2.683743386198e-04, -2.834525118783,  2.680023156066, -0.845229662945}
,{2.608434969930e-04, -2.831874866142,  2.674691791984, -0.842556082345}
,{2.204280700390e-04, -2.844713412111,  2.699061342074, -0.854127501894}
,{2.005004099728e-04, -2.849466946743,  2.707980937905, -0.858313490752}
,{1.827855384809e-04, -2.854053692629,  2.716612275857, -0.862375797689}
,{1.670003236655e-04, -2.858445294112,  2.724895751566, -0.866283457130}
,{1.528153586680e-04, -2.862665823265,  2.732873123563, -0.870054484939}
,{1.429801317530e-04, -2.865243505329,  2.737716241189, -0.872329755728}
,{1.323701711751e-04, -2.868502111494,  2.743876496279, -0.875242014614}
,{1.225311050362e-04, -2.871799123745,  2.750126888322, -0.878205233472}
,{1.114257076081e-04, -2.876231141130,  2.758575819195, -0.882233252358}
,{1.030159413657e-04, -2.879533439652,  2.764865228458, -0.885228772865}
,{9.508550482806e-05, -2.882867314570,  2.771226367634, -0.888263967559}
,{8.624390641154e-05, -2.887213699604,  2.779555660485, -0.892255716974}
,{8.230180925639e-05, -2.888410549470,  2.781807387295, -0.893314536016}
,{7.277858574106e-05, -2.893967410179,  2.792499242289, -0.898459053524}
,{7.119056581895e-05, -2.893881448938,  2.792283724840, -0.898331085336}
,{6.528133087447e-05, -2.897216513924,  2.798691601124, -0.901409805869}
,{6.316871578838e-05, -2.897808298469,  2.799796666492, -0.901925199307}
,{6.226901753736e-05, -2.896653005107,  2.797489168476, -0.900773894351}
,{5.545198290269e-05, -2.902312878640,  2.808453763165, -0.906085432542}
,{5.268715696216e-05, -2.903792467211,  2.811288832818, -0.907443678451}
,{4.870664319612e-05, -2.906584225062,  2.816675834737, -0.910042903032}
,{4.568271318928e-05, -2.908648965891,  2.820656388935, -0.911961740331}
,{4.363880773794e-05, -2.909891190474,  2.823043129832, -0.913108300551}
,{4.433166932649e-05, -2.907394571728,  2.818131140716, -0.910692237319}
,{3.861946150041e-05, -2.913647622573,  2.830297184907, -0.916610942872}
,{3.701831639158e-05, -2.914711907270,  2.832345375177, -0.917596449591}
,{3.759439057982e-05, -2.912349314702,  2.827692375329, -0.915305466236}
,{3.291917306925e-05, -2.918153615640,  2.839005212598, -0.920818677785}
,{3.132863504363e-05, -2.919475130046,  2.841560960072, -0.922054501391}
,{2.977356706657e-05, -2.920845804786,  2.844214822685, -0.923339244331}
,{2.852685306376e-05, -2.921880444782,  2.846214867792, -0.924305896157}
,{2.579470194231e-05, -2.925250194957,  2.852779963392, -0.927503973733}
,{2.502298736051e-05, -2.925586310490,  2.853414135397, -0.927802801919}
,{2.407727299691e-05, -2.926558555119,  2.855301375604, -0.928718743212}
,{2.304203791381e-05, -2.927592038189,  2.857305876026, -0.929690795799}
,{2.210365559085e-05, -2.928550780599,  2.859166021731, -0.930593137476}
,{2.098698274833e-05, -2.929866570638,  2.861725800548, -0.931838242927}
,{2.028452801639e-05, -2.930512719798,  2.862975927009, -0.932442922683}
,{2.116746309189e-05, -2.927694632010,  2.857428874535, -0.929713075062}
,{2.097961636860e-05, -2.927676206891,  2.857384635422, -0.929687448916}
,{1.918631005948e-05, -2.930143280582,  2.862194597064, -0.932032130173}
,{1.889364132801e-05, -2.930299399728,  2.862491165131, -0.932172871762}
,{1.809295032373e-05, -2.931314729577,  2.864467644761, -0.933134822234}
,{1.738388321115e-05, -2.932222168037,  2.866234197369, -0.933994645449}
,{1.669494877610e-05, -2.933134177538,  2.868010511213, -0.934859638726}
,{1.578385723255e-05, -2.934480827857,  2.870637855595, -0.936141243880}
,{1.256752862744e-05, -2.941856104391,  2.885113558721, -0.943244886801}
,{1.268309193192e-05, -2.941184963525,  2.883786915089, -0.942589268472}
,{1.287786374427e-05, -2.940379192085,  2.882197432927, -0.941805362978}
,{1.265662476801e-05, -2.940487078508,  2.882400270152, -0.941900535019}
,{1.187148436310e-05, -2.942054316368,  2.885470493159, -0.943404305307}
,{1.164247070573e-05, -2.942178952907,  2.885705881982, -0.943515286604}
,{1.118289582480e-05, -2.943047974034,  2.887406764328, -0.944347607399}
,{1.096996599548e-05, -2.943214692791,  2.887726394264, -0.944500731508}
,{1.051192676638e-05, -2.944098081132,  2.889455608752, -0.945347015693}
,{1.034431999514e-05, -2.944180238067,  2.889609645334, -0.945419062947}
,{8.899826335895e-06, -2.947949451173,  2.897014927052, -0.949056576053}
,{8.969457924612e-06, -2.947476340020,  2.896078936432, -0.948593626954}
,{8.913793207999e-06, -2.947370875232,  2.895866199282, -0.948486410256}
,{8.914095242840e-06, -2.947028073071,  2.895185269485, -0.948148282318}
,{9.223396251312e-06, -2.945244061251,  2.891659308833, -0.946406024186}
,{9.178148260069e-06, -2.945129312775,  2.891429099853, -0.946290608929}
,{8.612458924473e-06, -2.946530053737,  2.894178332352, -0.947639666156}
,{8.856724073802e-06, -2.945568196448,  2.892282574876, -0.946705521703}
,{8.226554968904e-06, -2.947055278548,  2.895199239180, -0.948135734077}
,{8.023233667331e-06, -2.947623554202,  2.896315883650, -0.948684306215}
,{8.064365332294e-06, -2.947239496147,  2.895555780544, -0.948308220032}
,{7.504412957715e-06, -2.948630347088,  2.898284996926, -0.949647145424}
,{7.128843179105e-06, -2.949810291044,  2.900606148641, -0.950788728754}
,{7.069271754645e-06, -2.949639322376,  2.900263004204, -0.950616612556}
,{6.707055331456e-06, -2.950827267180,  2.902600727285, -0.951766753049}
,{6.657745527883e-06, -2.950648170833,  2.902241944142, -0.951587115563}
,{6.322301762562e-06, -2.951791401786,  2.904492381156, -0.952694657068}
,{6.306992653249e-06, -2.951490682447,  2.903893995036, -0.952397005597}
,{5.980275953443e-06, -2.952639840745,  2.906156556435, -0.953510735414}
,{5.124202458906e-06, -2.956467151597,  2.913709100873, -0.957236825074}
,{5.038862653794e-06, -2.956625294648,  2.914017680606, -0.957387347095}
,{4.989129439248e-06, -2.956633506678,  2.914030529413, -0.957392033605}
,{5.011029182933e-06, -2.956313356249,  2.913395349174, -0.957076981896}
,{4.641222850266e-06, -2.957795911540,  2.916318014112, -0.958517461349}
,{4.564157959575e-06, -2.957952842757,  2.916624640991, -0.958667234076}
,{4.484693013951e-06, -2.958174316654,  2.917059486789, -0.958880685442}
,{4.394531719076e-06, -2.958365647726,  2.917433538931, -0.959063496673}
,{4.221147999806e-06, -2.959087077476,  2.918856087277, -0.959764788652}
,{4.077814021120e-06, -2.959638201001,  2.919941982459, -0.960299703644}
,{4.004439214444e-06, -2.959810988536,  2.920280246279, -0.960465253304}
,{3.870781167925e-06, -2.960347533710,  2.921337782606, -0.960986378115}
,{3.795852483002e-06, -2.960531787050,  2.921698666882, -0.961163083980}
,{3.657724093498e-06, -2.961105189567,  2.922829107827, -0.961720260537}
,{3.594593495970e-06, -2.961259247676,  2.923130756464, -0.961867914194}
,{3.470440250730e-06, -2.961811966731,  2.924221062164, -0.962405624992}
,{3.409227550266e-06, -2.961967584394,  2.924525902220, -0.962554908599}
,{3.348798692571e-06, -2.962122223280,  2.924828841377, -0.962703269298}
,{3.299539341528e-06, -2.962270597428,  2.925120217042, -0.962846320074}
,{3.195725366556e-06, -2.962614019682,  2.925795361896, -0.963178146489}
,{3.142288356184e-06, -2.962785503271,  2.926132351121, -0.963343705562}
,{3.108778021121e-06, -2.962889791000,  2.926337205906, -0.963444306128}
,{3.095672802989e-06, -2.962799821817,  2.926157068699, -0.963354151208}
,{2.918458531398e-06, -2.963718287838,  2.927970674588, -0.964249468291}
,{2.847490601310e-06, -2.964017238856,  2.928559857360, -0.964539771014}
,{2.868005573808e-06, -2.963829880317,  2.928188614442, -0.964355866119}
,{2.697173733024e-06, -2.964704402226,  2.929915136128, -0.965208036728}
,{2.724146780042e-06, -2.964447751225,  2.929406486777, -0.964956011405}
,{2.688792780425e-06, -2.964533241845,  2.929573671274, -0.965037740637}
,{2.408449694458e-06, -2.966309613891,  2.933086201742, -0.966774179401}
,{2.422470399277e-06, -2.966102603985,  2.932675198083, -0.966570171628}
,{2.422910418742e-06, -2.965936976614,  2.932345403833, -0.966406004309}
,{2.448045690273e-06, -2.965665860184,  2.931807790746, -0.966139482517}
,{2.356697305195e-06, -2.966202949739,  2.932869134750, -0.966663828314}
,{2.362469793682e-06, -2.966018621617,  2.932502611504, -0.966481627417}
,{2.266132158590e-06, -2.966621254536,  2.933693970847, -0.967070450179}
,{2.261575199336e-06, -2.966515281860,  2.933482532065, -0.966964988630}
,{2.168548806081e-06, -2.967117049321,  2.934672420598, -0.967553202727}
,{2.137025500226e-06, -2.967272979178,  2.934980061558, -0.967704945354}
,{2.124600705034e-06, -2.967209087022,  2.934851668713, -0.967640457090}
,{2.084914563993e-06, -2.967405764392,  2.935239681933, -0.967831832626}
,{2.031948978520e-06, -2.967707763446,  2.935836196869, -0.968126401474}
,{2.018330699238e-06, -2.967680213968,  2.935780060263, -0.968097827964}
,{1.971386640709e-06, -2.967965695600,  2.936344212515, -0.968376545528}
,{1.927300994464e-06, -2.968231690764,  2.936869812370, -0.968636194305}
,{1.881642727408e-06, -2.968511325225,  2.937422415054, -0.968909208186}
,{1.838171223234e-06, -2.968770010977,  2.937933487488, -0.969161638340}
,{1.784786005210e-06, -2.969126316449,  2.938638030987, -0.969509929752}
,{1.770777280918e-06, -2.969156900019,  2.938697583306, -0.969538912510}
,{1.721429112034e-06, -2.969483542193,  2.939343404495, -0.969858140873}
,{1.704506255806e-06, -2.969540685023,  2.939455561954, -0.969913172425}
,{1.629870453601e-06, -2.970172731312,  2.940707260279, -0.970532899097}
,{1.622623317044e-06, -2.970092835991,  2.940547340148, -0.970452881534}
,{1.581129134620e-06, -2.970415852809,  2.941186679549, -0.970769245611}
,{1.543582360597e-06, -2.970689168327,  2.941727389889, -0.971036677980}
,{1.490872660215e-06, -2.971079834249,  2.942500329729, -0.971419004608}
,{1.493007106945e-06, -2.970960527537,  2.942262824031, -0.971300803486}
,{1.436763278129e-06, -2.971429050857,  2.943190503178, -0.971760015558}
,{1.441720169471e-06, -2.971291258559,  2.942916454510, -0.971623754230}
,{1.438254784381e-06, -2.971227605566,  2.942789259052, -0.971560215232}
,{1.509614413919e-06, -2.970109814039,  2.940569477030, -0.970458153377}
,{1.480835536527e-06, -2.970304045938,  2.940953452780, -0.970647926006}
,{1.454572423842e-06, -2.970480412637,  2.941302091152, -0.970820223942}
,{1.438115132091e-06, -2.970573613032,  2.941486064463, -0.970911013316}
,{1.417372621937e-06, -2.970706746212,  2.941749137812, -0.971040974227}
,{1.399422786430e-06, -2.970817556527,  2.941968027551, -0.971149071602}
,{1.358983787392e-06, -2.971134851574,  2.942595901932, -0.971459691374}
,{1.332056797643e-06, -2.971287572706,  2.942897334989, -0.971608430226}
,{1.333964799399e-06, -2.971268863982,  2.942860265422, -0.971590067475}
,{1.277855308235e-06, -2.971748710796,  2.943810285959, -0.972060297307}
,{1.265111217896e-06, -2.971823816108,  2.943958563122, -0.972133481903}
,{1.242740730434e-06, -2.971993217396,  2.944293670342, -0.972299210206}
,{1.224752293805e-06, -2.972123149397,  2.944550607943, -0.972426233794}
,{1.214023388485e-06, -2.972183227280,  2.944669151035, -0.972484709732}
,{1.122120163788e-06, -2.973448090908,  2.947179723026, -0.973730509999}
,{1.147704452031e-06, -2.972849187692,  2.945988940143, -0.973138604747}
,{1.037141743376e-06, -2.974320672576,  2.948909141514, -0.974587431797}
,{1.026392722481e-06, -2.974386599589,  2.949039294550, -0.974651668568}
,{1.030055402418e-06, -2.974307633515,  2.948882320326, -0.974573656755}
,{1.025298522750e-06, -2.974276368975,  2.948819441934, -0.974542047661}
,{1.008420423498e-06, -2.974423464537,  2.949110663061, -0.974686190103}
,{9.918995742542e-07, -2.974566167670,  2.949393165441, -0.974826005872}
,{9.747305994567e-07, -2.974714518649,  2.949686843251, -0.974971349872}
,{9.590611657062e-07, -2.974855903556,  2.949966808496, -0.975109945879}
,{9.423116447138e-07, -2.975011678552,  2.950275324480, -0.975262703615}
,{9.329314341056e-07, -2.975058559987,  2.950367649002, -0.975308156083}
,{9.143289267932e-07, -2.975257154180,  2.950761289896, -0.975503221386}
,{9.109707178601e-07, -2.975174400036,  2.950595909975, -0.975420598968}
,{8.879370748671e-07, -2.975489349150,  2.951220955781, -0.975730718694}
,{8.793042625133e-07, -2.975543879173,  2.951328604044, -0.975783845567}
,{8.605752016555e-07, -2.975744653370,  2.951726555007, -0.975981041062}
,{8.432765128408e-07, -2.975931008636,  2.952095931988, -0.976164080076}
,{8.482585455161e-07, -2.975723382761,  2.951682615077, -0.975958384057}
,{8.220022250649e-07, -2.976116250797,  2.952462582600, -0.976345509801}
,{8.085320876150e-07, -2.976254499686,  2.952736515655, -0.976481207437}
,{7.944871556198e-07, -2.976406499847,  2.953037788757, -0.976630494424}
,{7.904375614043e-07, -2.976400756004,  2.953025819554, -0.976624273112}
,{7.780241284561e-07, -2.976527346674,  2.953276631116, -0.976748506418}
,{7.661672911663e-07, -2.976651442500,  2.953522535814, -0.976870327146}
,{7.542470839006e-07, -2.976778167776,  2.953773671395, -0.976994749372}
,{7.370211903002e-07, -2.976988859304,  2.954191536095, -0.977201939770}
,{7.340827722357e-07, -2.976979550599,  2.954172590358, -0.977192305676}
,{7.229753320459e-07, -2.977101713637,  2.954414724257, -0.977312287645}
,{7.117695792758e-07, -2.977225010094,  2.954659101760, -0.977433379896}
,{7.015101719521e-07, -2.977333691218,  2.954874457306, -0.977540064578}
,{6.871565074729e-07, -2.977515214068,  2.955234515938, -0.977718614713}
,{6.825414480272e-07, -2.977538951493,  2.955281232785, -0.977741598750}
,{6.718829236307e-07, -2.977660429421,  2.955522042926, -0.977860941622}
,{6.624451266424e-07, -2.977765968446,  2.955731229635, -0.977964598744}
,{6.525990766670e-07, -2.977881167809,  2.955959622419, -0.978077802011}
,{6.442779872939e-07, -2.977973335535,  2.956142289996, -0.978168310183}
,{6.349042187548e-07, -2.978081372934,  2.956356458499, -0.978274450661}
,{6.221827334452e-07, -2.978253555335,  2.956698085524, -0.978443908006}
,{6.196534713121e-07, -2.978247256957,  2.956685174462, -0.978437297851}
,{6.095745144208e-07, -2.978372130036,  2.956932811507, -0.978560071897}
,{5.972985356317e-07, -2.978546364202,  2.957278581896, -0.978731620396}
,{5.941410002874e-07, -2.978552288670,  2.957289954344, -0.978737071534}
,{5.859028012756e-07, -2.978657562833,  2.957498753134, -0.978840604398}
,{5.736663596378e-07, -2.978831405754,  2.957843738668, -0.979011759247}
,{5.729986893854e-07, -2.978797552138,  2.957776124472, -0.978977999335}
,{5.628244237332e-07, -2.978939117257,  2.958057020508, -0.979117340426}
,{5.538301909302e-07, -2.979064106245,  2.958305021514, -0.979240361439}
,{5.372390913250e-07, -2.979324454530,  2.958821902744, -0.979496910976}
,{5.306204149402e-07, -2.979411595573,  2.958994750922, -0.979582624728}
,{5.228446342542e-07, -2.979518377675,  2.959206604270, -0.979687703749}
,{5.168505597952e-07, -2.979590952529,  2.959350484944, -0.979759015565}
,{5.099638999972e-07, -2.979685428581,  2.959537917992, -0.979851979448}
,{5.042854714921e-07, -2.979768999472,  2.959703776310, -0.979934272553}
,{4.970399550830e-07, -2.979857534728,  2.959879300430, -0.980021268662}
,{4.913572847487e-07, -2.979948120481,  2.960059148006, -0.980110536168}
,{4.839647774757e-07, -2.980039466932,  2.960240252072, -0.980200301176}
,{4.784423931481e-07, -2.980129523047,  2.960419064118, -0.980289062628}
,{4.707752017108e-07, -2.980228570368,  2.960615482737, -0.980386441594}
,{4.647276138359e-07, -2.980311224130,  2.960779443166, -0.980467754308}
,{4.553137473717e-07, -2.980472101638,  2.961098933554, -0.980626376602}
,{4.571527634223e-07, -2.980414857764,  2.960985020529, -0.980569705612}
,{4.762188368090e-07, -2.979693001206,  2.959547879924, -0.979854402499}
,{4.784912976952e-07, -2.979556940363,  2.959276865533, -0.979719446679}
,{4.672062102751e-07, -2.979789613614,  2.959739237851, -0.979949157031}
,{4.617837122689e-07, -2.979903057301,  2.959964688361, -0.980061169276}
,{4.564684938924e-07, -2.980050980586,  2.960259002008, -0.980207564953}
,{4.499032248972e-07, -2.980146879660,  2.960449287509, -0.980301957946}
,{4.456890879689e-07, -2.980169343622,  2.960493425918, -0.980323636607}
,{4.400455917430e-07, -2.980283173683,  2.960719642797, -0.980436029068}
,{4.354376245130e-07, -2.980327419217,  2.960807175095, -0.980479320440}
,{4.289831673221e-07, -2.980458444556,  2.961067562746, -0.980608689206}
,{4.247494341802e-07, -2.980486748004,  2.961123340161, -0.980636167408}
,{4.196189097305e-07, -2.980589665845,  2.961327860673, -0.980737775209}
,{4.148948209437e-07, -2.980648283692,  2.961444043917, -0.980795345330}
,{4.099446451189e-07, -2.980734003570,  2.961614270837, -0.980879857322}
,{4.042743575461e-07, -2.980831752702,  2.961808376061, -0.980976219085}
,{4.000514436386e-07, -2.980881485900,  2.961906908428, -0.981025022477}
,{3.954297951303e-07, -2.980965794355,  2.962074369764, -0.981108179979}
,{3.904143468514e-07, -2.981049417522,  2.962240391944, -0.981190584008}
,{3.868626436976e-07, -2.981087708224,  2.962316200174, -0.981228105087}
,{3.816589664307e-07, -2.981200837698,  2.962541071047, -0.981339851690}
,{3.792820862936e-07, -2.981186962439,  2.962512973868, -0.981325632147}
,{3.782011330777e-07, -2.981174739879,  2.962488405860, -0.981313287780}
,{3.632580074520e-07, -2.981478805700,  2.963092580102, -0.981613411144}
,{3.681866533878e-07, -2.981349910946,  2.962836245036, -0.981485965903}
,{3.609823499007e-07, -2.981487813191,  2.963110184056, -0.981622009884}
,{3.562475395791e-07, -2.981572685519,  2.963278730990, -0.981705689224}
,{3.526817496935e-07, -2.981632379758,  2.963397238246, -0.981764505807}
,{3.491732788463e-07, -2.981691425812,  2.963514460744, -0.981822685759}
,{3.454738374264e-07, -2.981755557959,  2.963641797593, -0.981885894161}
,{3.439384185233e-07, -2.981770883106,  2.963672120171, -0.981900893127}
,{3.430944911509e-07, -2.981770631371,  2.963671480713, -0.981900506248}
,{3.363914724552e-07, -2.981909706437,  2.963947834726, -0.982037791898}
,{3.355022508700e-07, -2.981906814998,  2.963941915850, -0.982034765350}
,{3.285535872255e-07, -2.982053659594,  2.964233725652, -0.982179737505}
,{3.275112444578e-07, -2.982055780697,  2.964237781064, -0.982181672856}
,{3.209815632177e-07, -2.982201145008,  2.964526710425, -0.982325244436}
,{3.195858265803e-07, -2.982206733932,  2.964537617768, -0.982330564250}
,{3.139243140993e-07, -2.982348519265,  2.964819565071, -0.982470731881}
,{3.115608571136e-07, -2.982364040396,  2.964850108624, -0.982485756667}
,{3.110125621486e-07, -2.982358953025,  2.964839866648, -0.982480602610}
,{3.047177451476e-07, -2.982498742853,  2.965117700722, -0.982618653151}
,{3.039852910502e-07, -2.982495631065,  2.965111364737, -0.982615429687}
,{2.780038519390e-07, -2.983382545937,  2.966876668243, -0.983493844302}
};

static void
recgauss_3coeff_nonneg ( double sigma, double *a )
{
  int i = rint(sigma);
  if( i < 0 ) i = 0;
  if( i > RECGAUSS_MAX_SIGMA ) i = RECGAUSS_MAX_SIGMA;
  a[0] = recgauss_3coeff_nonneg_table[i][0];
  a[1] = recgauss_3coeff_nonneg_table[i][1];
  a[2] = recgauss_3coeff_nonneg_table[i][2];
  a[3] = recgauss_3coeff_nonneg_table[i][3];
}

// Find a quantile of weighted data
// 
// input:
// x[0],...,x[n-1] and corresponding non-negative weights w[0],...,w[n-1]
// The data is not sorted (otherwise the problem is trivial)
// 
// return x[k], such that sum_i^k w[i] <= p < sum_i^{k+1} w[i]
// The indices are such that the cumulative sums correspond to those
// of sorted data (although only partial quicksort is done here). 
// Both x and w will be scrambled by this routine.
// 
// example: p = sum_w/2 is the weighted median
//
// notes:
// - the weights don't have to sum to one; but if it does, p is the
//   cumulative probability
// - the function x(p), where p is the cumulative sum of weights,
//   is assumed to be a step function defined on intervals
//      [-inf,p_x0),[p_x0,p_x1),[p_x1,p_x2),...[p_x{n-1},+inf]
//   Thus, the larger x is returned if p is exactly on an interval boundary
// 
static double
wquantile ( 
    int n, double *x, // data
    double *w,        // weights corresponding to x, w[i] >= 0
    double sum_w,     // sum w[i]
    double p          // requested cumulative sum
    )
{
  double t, v;
  int a = 0, b = n-1;
  double sa = 0, sb = sum_w;
  while( b > a )
    {
    // put the midpoint of three at b
    if( x[b] < x[a] )
      {
      t = x[a]; x[a] = x[b]; x[b] = t; 
      v = w[a]; w[a] = w[b]; w[b] = v; 
      }
    int c = (a+b)/2;
    if( x[c] < x[b] )
      {
      t = x[c]; x[c] = x[b];
      v = w[c]; w[c] = w[b];
      if( t < x[a] )
        {
        x[b] = x[a]; x[a] = t;
        w[b] = w[a]; w[a] = v;
        }
      else
        { x[b] = t; w[b] = v; }
      }
    
    int i = a-1, j = b;
    double sl = sa, sr = sb;
    for(;;)
      {
      while( x[++i] < x[b] )
        sl += w[i];
      while( x[b] < x[--j] )
        sr -= w[j];
      if( i >= j ) break;
      t = x[i]; x[i] = x[j]; x[j] = t;
      v = w[i]; w[i] = w[j]; w[j] = v;
      sl += w[i]; sr -= w[j];
      }
    t = x[i]; x[i] = x[b]; x[b] = t;
    v = w[i]; w[i] = w[b]; w[b] = v;
    
    if( p < sl )
      { b = i - 1; sb = sl; }
    else if( p >= sr )
      { a = i + 1; sa = sr; }
    else
      return x[i];
    }
  return x[b];
}


/* added by Mark Robinson, 23rd March 2008 to interface with R */

void call_loexp(int *intparams, double *dblparams, double *y, double *w, double *Ey, double *ow) {
   int ret=0,i;
   ret = loexp( intparams[1], y, w, intparams[2], dblparams[0], dblparams[1], dblparams[2], dblparams[3], intparams[3], Ey, ow);
   intparams[0] = ret;
}


// 
// loexp() fits local expectiles.
//
// The input is series of (x0,y0),(x1,y1),...  However, x_i - x_{i-1}
// is assumed constant, and the x values are not considered. Each observation
// can be given non-negative weight. This routine find 
//   E(y) = beta0 + beta1 x + beta2 x^2
// at fitted locally each point, using Gaussian kernel to weigh nearby
// points. The kernel width is determined by sigma (defined in the number
// of data points).
//
// beta's for each data point i were found to minimize
//
//   \sum_{j=1}^n w_j G_{ij} R_j A_j (y_j-E[y|x_j,beta_i])^2
// 
// which is local (G is gaussian kernel), robust (R is Tukey's biweight 
// determined by the median of residuals) and asymmetric (A is the expectile
// weight: 1-alpha if the residual is negative and alpha if positive).
// 
// The robustness and asymmetry weights are based on pooled residuals
// across all data points (they are not local for a given i). 
// However, the robustness weights are determined separately for
// positive and negative residuals.
// 
// The polynomial dummy x's is arbitrary and thus beta's are not
// meaningful and not returned. Ey is invariant to the choice of dummies.
//
// return:
//  0  converge
//  1  iteration limit exceeded
//  2  non-positive definite X'WX encoutered when solving local regression
//  3  invalid order of polynomials
//  
int
loexp (
    int n,                // data length
    const double *y,      // data vector
    const double *w,      // sample weights (assume w[i]=1 if w == NULL)
    int polynomial_order, // 0, 1 or 2
    double sigma,         // gaussian kernel SD (in # of data points)
    double alpha,         // asymmetry weight [0..1]
    double biweight,      // Tukey's biweight tuning

    double tolerance,     // convergence tolerance
    int iter_max,         // maximum iteration
    
    // output (arrays should be allocated by caller)
    double *Ey,           // fitted curve
    double *v_            // final weights (w_biweight * w_asymmetry),
                          // can be NULL if not needed
    )
{

  /*
  printf("------------------------\nParameters:\n");
  printf("length=%d\n",n);
  printf("polynomial_order=%d\n",polynomial_order);
  printf("iter_max=%d\n",iter_max);
  printf("sigma=%f\n",sigma);
  printf("alpha=%f\n",alpha);
  printf("biweight=%f\n",biweight);
  printf("tolerance=%f\n",tolerance);
  printf("------------------------\n");
  */


  const int RIGHT_PAD = 8; // 8*sigma points needed to avoid kernel distortion
  
  const int K = 3;                       // filter coefficients
  int nn = K + n + RIGHT_PAD * sigma;
  
  double *workspace = alloc(nn*8 + 3*n); // 8 is # sums in quadratic regression
  double *v = workspace + nn*8;          // working regression weights
  double *AR = v + n;                    // absolute residual
  double *wAR = AR + n;                  // weights of absolute residual
  set_d ( 8*K, 0, workspace );           // unused left padding

  if( alpha < 0 ) alpha = 0;
  if( alpha > 1 ) alpha = 1;

  double filter_coeff[4];
  recgauss_3coeff_nonneg ( sigma, filter_coeff );  

  set_d ( n, 1, v ); // start with OLS
  
  /* int i;
  for(i=0; i<n; i++) {
    printf("v[%d]=%f\n",i,v[i]);
  } */
  
  double RSS = 0, oldRSS = 0.0;
  int return_status = 1;
  for(int r = 0; r < iter_max; r++ )
    {
	/* printf("Iteration %d.\n",r); */
    // fill in the LHS-RHS of local regression
    double *s = workspace + 8*K;
    for(int i = 0; i < n; i++, s += 8 )
      {
      double x = 4.0/n * i - 2.0; // polynomial dummy vars
      s[0] = v[i];
      s[1] = v[i]*x;
      s[2] = s[1]*x;
      s[3] = s[2]*x;
      s[4] = s[3]*x;
      s[5] = v[i]*y[i];
      s[6] = s[5]*x;
      s[7] = s[6]*x;
      }
    set_d ( RIGHT_PAD * sigma * 8, 0.0, s );
	
	/* int i; 
    for(i=0; i<n*8; i++) {
      printf("s[%d]=%f\n",i,s[i]);
    } 

	printf("After set_d.\n"); */
    
    // smooth them
    recgauss_filter_col( nn, 8, workspace, 3, filter_coeff );
	//printf("After recgauss_filter_col (polynomial_order=%d).\n",polynomial_order);

    // solve and find Ey
    s = workspace + 8*K;
    RSS = 0;
    int ineg = 0, ipos = n; // index to negative and positive residuals
    double swneg = 0, swpos = 0;
    for(int i = 0; i < n; i++, s += 8 )
      {
	  //printf("%d.\n",i);
      double sw = s[0], sx1 = s[1], sx2 = s[2], sx3 = s[3], sx4 = s[4],
             sy = s[5], sxy = s[6], sx2y = s[7];
      double beta0 = 0, beta1 = 0, beta2 = 0;
      if(polynomial_order == 2 )
        {
        double det = sx2*(2*sx1*sx3 + sw*sx4 - sx2*sx2)
                     - sw*sx3*sx3 -sx1*sx1*sx4;
        if( det <= 0 )
          { return_status = 2; goto RETURN; }
        beta0 = ( sx2*(sx3*sxy + sx4*sy - sx2*sx2y)
                  + sx1*(sx3*sx2y - sx4*sxy) - sx3*sx3*sy )/det;
        beta1 = ( sx2*(sx1*sx2y + sx3*sy - sx2*sxy)
                  + sx4*(sw*sxy - sx1*sy) - sw*sx3*sx2y )/det;
        beta2 = ( sx2*(sw*sx2y + sx1*sxy - sx2*sy)
                  + sx1*(sx3*sy - sx1*sx2y) -sw*sx3*sxy )/det;
        }
      else if(polynomial_order == 1 )
        {
        double det = sw*sx2 - sx1*sx1; 
        if( det <= 0 )
          { return_status = 2; goto RETURN; }
        beta0 = (sx2*sy - sx1*sxy)/det;
        beta1 = (sw*sxy - sx1*sy)/det;
        beta2 = 0;
        }
      else if( polynomial_order == 0 )
        {
        if( sw <= 0 )
          { return_status = 2; goto RETURN; }
        beta0 = sy/sw; beta1 = beta2 = 0;
        }
      else
        { return_status = 3; goto RETURN; }

      double x = 4.0/n * i - 2.0; // polynomial dummy vars
      Ey[i] = beta0 + beta1*x + beta2*x*x;
      double res = y[i] - Ey[i];
      RSS += v[i] * res * res;
      v[i] = res;
      
      double wres = w ? w[i] : 1;
      if( wres <= 0 ) continue;
      if( res < 0 )
        {
        swneg += wres;
        wAR[ineg] = wres;
        AR[ineg++] = fabs(res);
        }
      else
        {
        swpos += wres;
        AR[--ipos] = res;
        wAR[ipos] = wres;
        }
      }
/*
#ifdef LOEXP_SHELL
    fprintf(stderr,"%d RSS = %g\n", r, RSS);
#endif
*/
    if( r > 0 && 2*fabs(RSS - oldRSS)/fabs(RSS+oldRSS) < tolerance )
      { return_status = 0; break; }
    oldRSS = RSS;

    // update the weight
    
    // the scale of robustness kernel is MAR/0.6745 
    // (ref: John Fox 2002, appendix to "An R and S-PLUS Companion to
    // Applied Regression"). 
    // The recommended biweight cutoff is 4.685 of the scale.
    // 
    double bwcutneg = biweight
                       * wquantile( ineg, AR, wAR, swneg, swneg/2 )/0.6745;
    double bwcutpos = biweight
           * wquantile( n-ipos, AR + ipos, wAR + ipos, swpos, swpos/2 )/0.6745;
    
    for(int i = 0; i < n; i++ )
      {
      double res = v[i];
      if( res < 0 )
        {
        if( -res < bwcutneg )
          {
          double u = res/bwcutneg;
          u = 1 - u*u;
          v[i] = (1-alpha) * u * u;
          }
        else
          v[i] = 0;
        }
      else
        {
        if( res < bwcutpos )
          {
          double u = res/bwcutpos;
          u = 1 - u*u;
          v[i] = alpha * u * u;
          }
        else
          v[i] = 0;
        }
      if( w ) v[i] *= w[i];
      if( v_ ) v_[i] = v[i];
      }
    }

RETURN:
  free(workspace);
  return return_status;
}


/* I've commented out  everything below this line  */


/*

#include <stdio.h>
#include <getopt.h>

#ifdef LOEXP_SHELL

static int
get_double ( char *filename, double **y_ )
{
  double *y = 0; 
  int nallocated = 0;
  int n = 0;
  FILE *F = 0;
  if( 0 == strcmp(filename,"-"))
    F = stdin;
  else
    F = fopen(filename,"r");
  if(!F) return 0;
  
  while(!feof(F))
    {
    if( n >= nallocated )
      y = realloc( y, sizeof(double)*(nallocated += 100));
    if( 1 != fscanf(F, " %lf", &y[n]) ) break;
    n++;
    }
  *y_ = y;
  fclose(F);
  return n;
}

static int
get_double2 ( char *filename, double **x_, double **y_ )
{
  double *x = 0, *y = 0; 
  int nallocated = 0;
  int n = 0;
  FILE *F = 0;
  if( 0 == strcmp(filename,"-"))
    F = stdin;
  else
    F = fopen(filename,"r");
  if(!F) return 0;
  
  while(!feof(F))
    {
    if( n >= nallocated )
      {
      y = realloc( y, sizeof(double)*(nallocated += 100));
      x = realloc( x, sizeof(double)*nallocated);
      }
    if( 2 != fscanf(F, " %lf %lf", &x[n], &y[n]) ) break;
    n++;
    }
  *x_ = x;
  *y_ = y;
  fclose(F);
  return n;
}

static void
usage(char *argv0)
{
  fprintf(stderr,
  "usage:\n"
  " %s [options] [ <input> [ <output> ] ]\n"
  "\n"
  "Input and output uses stdin and stdout if not specified\n"
  "\n"
  "Options:\n"
  "-h        this help\n"
  "-w        input has sample weights in the second column [no weights]\n"
  "-a <num>  asymmetry [0.5]\n"
  "-s <num>  sigma [40]\n"
  "-g <num>  Tukey's biweight robustness [4.685]\n"
  "-r <int>  order of polynomials (0,1, or 2), [2]\n"
  "-t <num>  tolerance [1e-4]\n"
  "-R <int>  maximum number of iterations [100]\n"
  "-v <file> output the product of robustness and asymmetry weights\n"
  
  "\n",
  argv0 );
  exit(1);
}


int
main(int argc, char* argv[])
{
  char *infile = "-";
  char *outfile = "-";
  int weighted = 0;
  double sigma = 40;
  double alpha = 0.5;
  char *output_weight = 0;
  int poly_order = 2;
  double tolerance = 1e-4;
  int iter_max = 100;
  double biweight = 4.685; // recommended cutoff

  for(;;)
    {
    static struct option option_def[] = 
      {
        {"help", 0, 0, 'h'},
        { 0, 0, 0, 0 }
      };
    int c = getopt_long ( argc, argv, "hws:a:g:v:r:t:R:", option_def, 0 );
    if( c == -1 ) break;
    switch(c)
      {
      case 'g':
        biweight = atof(optarg);
        break;
      case 't':
        tolerance = atof(optarg);
        break;
      case 'R':
        iter_max = atoi(optarg);
        break;
      case 'r':
        poly_order = atoi(optarg);
        if(poly_order < 0 || poly_order > 2 )
          usage(argv[0]);
        break;
      case 'v':
        output_weight = optarg;
        break;
      case 'a':
        alpha = atof(optarg);
        break;
      case 's':
        sigma = atof(optarg);
        break;
      case 'w':
        weighted = 1;
        break;
      case 'h':
      default:
        usage(argv[0]);
        break;
      }
    }
  if(argc - optind > 0) infile = argv[optind];
  if(argc - optind > 1) outfile = argv[optind+1];
   
  double *y = 0, *w = 0;
  int n;
  if( weighted )
    n = get_double2(infile, &y, &w );
  else
    n = get_double(infile, &y);

  if( n == 0 )
    { fprintf(stderr,"no data.\n"); exit(1); }
  else
    fprintf(stderr,"\ndata: n = %d %s\n",
        n, w ? "with weights":"without weights");

  double *v = alloc(n), *Ey = alloc(n);

  int ret = loexp( n, y, w, poly_order, sigma, alpha, 
      biweight, tolerance, iter_max,
      Ey, v);
  if( ret != 0 )
    fprintf(stderr,"loexp failed with code %d\n", ret);
  
  if(strcmp("-",outfile) )
    freopen(outfile,"w",stdout);
  for(int i = 0; i < n; i++ )
    fprintf(stdout, "%g\n",Ey[i]);

  if(output_weight)
    {
    freopen(output_weight,"w",stdout);
    for(int i = 0; i < n; i++ )
      fprintf(stdout, "%g\n",v[i]);
    }
  
  exit(0);
}
#endif

*/
