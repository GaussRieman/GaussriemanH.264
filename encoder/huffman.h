
#ifndef X264_0_148_R2795_1_HUFFMAN_H
#define X264_0_148_R2795_1_HUFFMAN_H

#include "common/common.h"
#define safe_free(p) do{if(p!=NULL){free(p);p=NULL;}}while(0)
typedef struct HuffmanCode {
    uint32_t value;
    uint32_t bit_length;
} HuffmanCode;

typedef struct INode {
    int left;
    int right;
    int value;
    int frequency;
    int index;
} INode;

typedef struct HuffmanTree
{
    int *pos;
    int root_pos;
    struct INode *iNode;
    int size;
}HuffmanTree;

typedef struct HuffmanEncodeContext {
    HuffmanCode *huffman_code;
    HuffmanTree *huffman_tree;
} HuffmanEncodeContext;

typedef struct HuffmanDecodeContext {
    HuffmanTree *huffman_tree;
} HuffmanDecodeContext;

typedef enum HuffCodeArrayType {
    CharArray, ShortArray, IntArray
} HuffCodeArrayType;



int initializeHuffmanTree(const unsigned int max_symbol_category,HuffmanTree **huffman_tree);

void releaseHuffmanTree(HuffmanTree *huffman_tree);


int initializeHuffmanEncodeContext(const unsigned int range, HuffmanEncodeContext **huffman_encode_context);

void releaseHuffmanEncodeContext(HuffmanEncodeContext *huffman_encode_context);

int encodeWithHuffman(const short *data, const int data_length, const int max_symbol_category,
                      HuffmanEncodeContext *huffman_encode_context, bs_t *bitstream);

#endif
