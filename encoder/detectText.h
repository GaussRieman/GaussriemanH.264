//
// Created by yanz on 2017/9/5.
//

#ifndef X264_0_148_R2795_1_DETECTTEXT_H
#define X264_0_148_R2795_1_DETECTTEXT_H
#include "common\common.h"

#define NGH_DIFF                                                                                      18
#define NBC_DIFF                                                                                       4
#define NGH_TH                                                                                         1
#define NBC_TH0                                                                                      179
#define NBC_TH1                                                                                      160

int nhgGenerateSmart_Ex(pixel *data, int stride, int length, int max_value);

int nbcGenerateSmart_Ex(pixel*data, int stride, int length, int max);

int decideTextBlock(pixel *p_dst, int stride_p_dst, pixel* p_src , int stride_p_src ,int length);

#endif //X264_0_148_R2795_1_DETECTTEXT_H
