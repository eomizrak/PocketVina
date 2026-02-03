/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "coords.h"
/*
	input¡ª¡ªa¡¢b£ºvec½á¹¹ÌåÏòÁ¿£»
	inter¡ª¡ªacc:double;
	demo:Ç°Ìáa¡¢bÔªËØ¸öÊýÒ»Ñù£¬²¢ÇÒ²»Îª¿Õ£¬·ñÔò·µ»Ø0
	a{£¨7£¬8£¬9£©£¬£¨5£¬6£¬7£©}£¬b{£¨1£¬2£¬3£©£¬£¨3£¬4£¬5£©} size=2
	acc=£¨(7-1)^2+(8-2)^2+(9-3)^2£©+£¨(5-3)^2+(6-4)^2+(7-5)^2£©=120
	·µ»Ø £¨120/2£©^(1/2)=7.74597
*/
fl rmsd_upper_bound(const vecv& a, const vecv& b) {
	VINA_CHECK(a.size() == b.size());
	fl acc = 0;
	VINA_FOR_IN(i, a)       //for(i=0;i<a.size;i++)
		acc += vec_distance_sqr(a[i], b[i]);
	return (a.size() > 0) ? std::sqrt(acc / a.size()) : 0;
}
fl rmsd_upper_bound_sqr(const vecv& a, const vecv& b) {
	VINA_CHECK(a.size() == b.size());
	fl acc = 0;
	VINA_FOR_IN(i, a)
		acc += vec_distance_sqr(a[i], b[i]);
	return (a.size() > 0) ? (acc / a.size()) : 0;
}



/*
	input¡ª¡ªa¡¢b.coords£ºvec½á¹¹ÌåÏòÁ¿£¬b£ºoutput_type½á¹¹ÌåÏòÁ¿
	output¡ª¡ªtmp£º<unint,1.79769e+308>,¿ÉÒÔÓÃfirst¡¢secondµ÷ÓÃ
	demo:Ç°Ìáa¡¢b.coordsÔªËØ¸öÊýÒ»Ñù£¬²¢ÇÒ²»Îª¿Õ£¬·ñÔò·µ»Øres=0
	 a{£¨7£¬8£¬9£©£¬£¨5£¬6£¬7£©}£¬b.coords{£¨1£¬2£¬3£©£¬£¨3£¬4£¬5£©} size=2
	res=[(£¨(7-1)^2+(8-2)^2+(9-3)^2£©+£¨(5-3)^2+(6-4)^2+(7-5)^2£©)/2]^(1/2)=7.74597
	·µ»Øtmp<2,7.74597>
*/

std::pair<sz, fl> find_closest(const vecv& a, const output_container& b) {
	std::pair<sz, fl> tmp(b.size(), max_fl);          //max_fl=1.79769e+308
	fl best_rmsd_sqr = max_fl;
	VINA_FOR_IN(i, b) { 
		fl res_sqr = rmsd_upper_bound_sqr(a, b[i].coords);
		if(i == 0 || res_sqr < best_rmsd_sqr) {
			best_rmsd_sqr = res_sqr;
			tmp.first = i;
		}
	}
	if(tmp.first < b.size() && best_rmsd_sqr < max_fl)
		tmp.second = std::sqrt(best_rmsd_sqr);
	return tmp;
}


/*
	input¡ª¡ªout£ºoutput_type½á¹¹ÌåÏòÁ¿£¬t£ºoutput_type½á¹¹Ìå£¬min_rmsd£ºdouble£¬max_size£ºunint
	¸ù¾Ýfind_closestº¯ÊýµÄdemo£¬closest_rmsd<2,7.74597>
	¡¤Èç¹ûoutÏòÁ¿ÔªËØ>2£¬min_rmsd>7.74597£¨±íÊ¾ÓÐÒ»¸ö·Ç³£Ð¡µÄÖµ£©²¢ÇÒt.e < out[2].e
		out[2]=t;
	¡¤·ñÔò£¨±íÊ¾Ã»ÓÐÐ¡µÄÖµ£©
	Èç¹ûoutÏòÁ¿ÔªËØ<max_size,ÔÚoutÏòÁ¿Ä©Î²¼ÓÉÏt£¨output_type£©½á¹¹Ìå
	·ñÔòÈç¹ûoutÏòÁ¿²»Îª¿Õ²¢ÇÒt.e < outÄ©Î²ÔªËØµÄeÖµ£¬½«outÄ©Î²ÔªËØ¸³ÖµÎªt
	¡¤×îºó½«outÏòÁ¿ÒÔout½á¹¹ÌåÖÐe½øÐÐÉýÐòÅÅÐò

*/
void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd, sz max_size) {
	if(max_size == 0) return;
	std::pair<sz, fl> closest_rmsd = find_closest(t.coords, out);
	if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
		if(t.e < out[closest_rmsd.first].e) { // the new one is better, apparently
			out[closest_rmsd.first] = t; // FIXME? slow
		}
	}
	else { // nothing similar
		if(out.size() < max_size)   //newå¼€è¾Ÿçš„ç©ºé—´åœ¨å †ä¸Šï¼Œè€Œä¸€èˆ¬å£°æ˜Žçš„å˜é‡å­˜æ”¾åœ¨æ ˆä¸Š
			out.push_back(new output_type(t));
		else if(!out.empty()) {
			sz worst_i = 0;
			fl worst_e = out[0].e;
			VINA_FOR_IN(i, out) {
				if(out[i].e > worst_e) {
					worst_e = out[i].e;
					worst_i = i;
				}
			}
			if(t.e < worst_e)
				out[worst_i] = t; // FIXME? slow
		}
	}
}
