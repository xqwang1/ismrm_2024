/* Copyright 2013-2014. The Regents of the University of California.
 * Copyright 2016-2017. Martin Uecker.
 * All rights reserved. Use of this source code is governed by
 * a BSD-style license which can be found in the LICENSE file.
 */

#include "misc/types.h"

typedef struct iter3_conf_s { TYPEID* TYPEID; } iter3_conf;

struct iter_op_s;



struct iter3_irgnm_conf {

	INTERFACE(iter3_conf);

	int iter;
	float alpha;
	float alpha_min;
	float alpha_min0;
	float redu;

	int cgiter;
	float cgtol;

	_Bool nlinv_legacy;
};





struct iter3_landweber_conf {

	INTERFACE(iter3_conf);

	int iter;
	float alpha;
	float epsilon;
};




extern const struct iter3_irgnm_conf iter3_irgnm_defaults;
extern const struct iter3_landweber_conf iter3_landweber_defaults;

