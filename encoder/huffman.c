#include "huffman.h"
#define HUFF_LEFT_VALUE 0  
#define HUFF_RIGHT_VALUE 1 

#define ByteLength 8
#define LEFT(x) (2 * (x) + 1)
#define RIGHT(x) (2 * (x) + 2)
#define PARENT(x) ((x) / 2)

int initializeHuffmanTree(const unsigned int max_symbol_category,HuffmanTree **huffman_tree){
    *huffman_tree=(HuffmanTree *) calloc(1, sizeof(HuffmanTree));
    if (NULL==(*huffman_tree)){
        return -1;
    }
    (*huffman_tree)->iNode=(INode *) calloc(max_symbol_category+max_symbol_category, sizeof(INode));
    if (NULL==(*huffman_tree)->iNode){
        releaseHuffmanTree(*huffman_tree);
        (*huffman_tree)=NULL;
        return -2;
    }
    (*huffman_tree)->pos=(int *)calloc(max_symbol_category, sizeof(int));
    if (NULL==(*huffman_tree)->pos){
        releaseHuffmanTree(*huffman_tree);
        (*huffman_tree)=NULL;
        return -3;
    }
    return 1;
}


void releaseHuffmanTree(HuffmanTree *huffman_tree){
    if (NULL != huffman_tree){
        safe_free(huffman_tree->iNode);
        safe_free(huffman_tree->pos);
        safe_free(huffman_tree);
    }
}

int getHuffmanCodeDict(HuffmanEncodeContext *huffman_encode_context, const int *frequencies, const unsigned int frequency_count);

static void codeHuffmanBitStream(const void *array, const int num_element,const HuffCodeArrayType array_element_type,
                                 const HuffmanCode *huffman_code_map, bs_t *out_stream);

static void codeHuffmanDictionary(struct HuffmanTree *huffman_tree, bs_t *out_stream, const int leaf_node_code_length);

static void buildTree(HuffmanTree *huffman_tree,const int *frequencies, const unsigned int frequency_count);

static void generateCodes(struct HuffmanTree *Huffman_tree,size_t node,const HuffmanCode *prefix, HuffmanCode *huffman_code_dict);

static inline void huffCodeAppendBit(int bit, HuffmanCode *huff_code);

static void codeHuffmanNodeBinary(struct HuffmanTree *huffman_tree, size_t node, bs_t *bit_stream, const int leaf_node_code_length);

static void generateShortFrequencyMap(const short *data, const int data_length, const int max_symbol_category, int *frequency);

static void addOneNode(struct HuffmanTree *huffman_tree,size_t index);

static INode popOneNode(struct HuffmanTree *huffman_tree);

static void adjustTreeHeapify(struct HuffmanTree *Huffman_tree,size_t idx);

int initializeHuffmanEncodeContext(const unsigned int max_encode_context, HuffmanEncodeContext **huffman_encode_context){
    *huffman_encode_context=(HuffmanEncodeContext *) calloc(1, sizeof(HuffmanEncodeContext));
    if (NULL==*huffman_encode_context){
        return -1;
    }
    (*huffman_encode_context)->huffman_code=(HuffmanCode *) calloc(max_encode_context, sizeof(HuffmanCode));
    if (NULL==(*huffman_encode_context)->huffman_code){
        releaseHuffmanEncodeContext(*huffman_encode_context);
        return -1;
    }
    HuffmanTree *huffman_tree=NULL;
    if (-1==initializeHuffmanTree(max_encode_context,&huffman_tree)){
        releaseHuffmanEncodeContext(*huffman_encode_context);
        return -1;
    }
    (*huffman_encode_context)->huffman_tree=huffman_tree;
    return 1;
}
void releaseHuffmanEncodeContext(HuffmanEncodeContext *huff_encode_context){
    if (NULL != huff_encode_context){
        safe_free(huff_encode_context->huffman_code);
        releaseHuffmanTree(huff_encode_context->huffman_tree);
        huff_encode_context->huffman_tree=NULL;
        safe_free(huff_encode_context);
    }
}

/**
* Adds a new element to the huffman tree .
*/
static void addOneNode(struct HuffmanTree *huffman_tree,size_t index){
    size_t i=huffman_tree->size;
    int tmp;
    huffman_tree->pos[huffman_tree->size]=index;
    huffman_tree->size++;
    while (i>0 && huffman_tree->iNode[huffman_tree->pos[i]].frequency<huffman_tree->iNode[huffman_tree->pos[PARENT(i)]].frequency){
        tmp=huffman_tree->pos[i];
        huffman_tree->pos[i]=huffman_tree->pos[PARENT(i)];
        huffman_tree->pos[PARENT(i)]=tmp;
        i = PARENT(i);
    }
}

/**
* Returns the element with the biggest priority from the tree .
*/
static INode popOneNode(struct HuffmanTree *huffman_tree){
    struct INode data;
    if (huffman_tree->size<1){
        /* Priority tree is empty */
        data.index=-1;
        return data;
    }
    data=huffman_tree->iNode[huffman_tree->pos[0]];
    huffman_tree->pos[0]=huffman_tree->pos[--huffman_tree->size];
    adjustTreeHeapify(huffman_tree,0);
    return data;
}

/**
* Turn an "almost-heap" into a heap .
*/
static void adjustTreeHeapify(struct HuffmanTree *Huffman_tree,size_t idx){
    size_t l_idx,r_idx,lrg_idx,tmp;
    l_idx = LEFT(idx);
    r_idx = RIGHT(idx);
    if (l_idx<Huffman_tree->size && Huffman_tree->iNode[Huffman_tree->pos[idx]].frequency>Huffman_tree->iNode[Huffman_tree->pos[l_idx]].frequency){
        lrg_idx=l_idx;
    }else{
        lrg_idx=idx;
    }
    if (r_idx<Huffman_tree->size && Huffman_tree->iNode[Huffman_tree->pos[r_idx]].frequency<Huffman_tree->iNode[Huffman_tree->pos[lrg_idx]].frequency){
        lrg_idx=r_idx;
    }
    if (lrg_idx!=idx){
        tmp=Huffman_tree->pos[idx];
        Huffman_tree->pos[idx]=Huffman_tree->pos[lrg_idx];
        Huffman_tree->pos[lrg_idx]=tmp;
        adjustTreeHeapify(Huffman_tree,lrg_idx);
    }
}
static void buildTree(HuffmanTree *huffman_tree,const int *frequencies, const unsigned int frequency_count) {
    huffman_tree->root_pos=0;
    huffman_tree->size=0;
    for (int i = 0; i < frequency_count; ++i) {
        if (frequencies[i] != 0) {
            huffman_tree->iNode[huffman_tree->size].value=i;
            huffman_tree->iNode[huffman_tree->size].frequency=frequencies[i];
            huffman_tree->iNode[huffman_tree->size].left=-1;
            huffman_tree->iNode[huffman_tree->size].right=-1;
            huffman_tree->iNode[huffman_tree->size].index=huffman_tree->size;
            addOneNode(huffman_tree,huffman_tree->size);
        }
    }
    size_t total_node=huffman_tree->size;

    if (total_node<1){
        printf("no huffman tree.\n");
        huffman_tree->root_pos=-1;
        return;
    }
    if (total_node==1){
        huffman_tree->iNode[huffman_tree->size].value=(frequencies[0])==0?0:1;
        huffman_tree->iNode[huffman_tree->size].frequency=0;
        huffman_tree->iNode[huffman_tree->size].left=-1;
        huffman_tree->iNode[huffman_tree->size].right=-1;
        huffman_tree->iNode[huffman_tree->size].index=huffman_tree->size;
        addOneNode(huffman_tree,huffman_tree->size);
        total_node=huffman_tree->size;
    }
    struct INode child_right,child_left;
    while (huffman_tree->size>1){
        child_right=popOneNode(huffman_tree);
        child_left=popOneNode(huffman_tree);
        huffman_tree->iNode[total_node].frequency=child_right.frequency+child_left.frequency;
        huffman_tree->iNode[total_node].left=child_left.index;
        huffman_tree->iNode[total_node].right=child_right.index;
        huffman_tree->iNode[total_node].index=total_node;
        addOneNode(huffman_tree,total_node);
        total_node++;
    }
    if (huffman_tree->size==1){
        huffman_tree->root_pos=total_node-1;
    }else{
        huffman_tree->root_pos=-1;
    }
}

static void generateCodes(struct HuffmanTree *Huffman_tree,size_t node,const HuffmanCode *prefix, HuffmanCode *huffman_code_dict) {
    if (Huffman_tree->iNode[node].left==-1 && Huffman_tree->iNode[node].right==-1){
        huffman_code_dict[Huffman_tree->iNode[node].value] = *prefix;
    }else{
        if (Huffman_tree->iNode[node].left!=-1){
            HuffmanCode left_prefix = *prefix;
            huffCodeAppendBit(HUFF_LEFT_VALUE, &left_prefix);
            generateCodes(Huffman_tree,Huffman_tree->iNode[node].left, &left_prefix, huffman_code_dict);
        }

        if (Huffman_tree->iNode[node].right!=-1) {
            HuffmanCode right_prefix = *prefix;
            huffCodeAppendBit(HUFF_RIGHT_VALUE, &right_prefix);
            generateCodes(Huffman_tree,Huffman_tree->iNode[node].right, &right_prefix, huffman_code_dict);
        }
    }
}


static inline void huffCodeAppendBit(int bit, HuffmanCode *huff_code) {
    huff_code->value = huff_code->value << 1 | bit;
    huff_code->bit_length += 1;
}

int getHuffmanCodeDict(HuffmanEncodeContext *huffman_encode_context, const int *frequencies, const unsigned int frequency_count) {

    buildTree(huffman_encode_context->huffman_tree,frequencies, frequency_count);
    if (huffman_encode_context->huffman_tree->root_pos!=-1){
        HuffmanCode prefix;
        prefix.bit_length = 0;
        prefix.value = 0;
        generateCodes(huffman_encode_context->huffman_tree,huffman_encode_context->huffman_tree->root_pos, &prefix, huffman_encode_context->huffman_code);
        return 1;
    }
    else{
        printf("Nil root detected.\n");
        return -1;
    }
}

static void codeHuffmanBitStream(const void *array, const int num_element, const HuffCodeArrayType array_element_type,
                                 const HuffmanCode *huffman_code_map,bs_t *bit_stream) {

    if (array_element_type == ShortArray) {
        for (int i = 0; i < num_element; i++) {
            HuffmanCode huffman_code = huffman_code_map[((short *) array)[i]];
            bs_write(bit_stream,huffman_code.bit_length,(uint32_t)huffman_code.value);
        }
    } else if (array_element_type == CharArray) {
        for (int i = 0; i < num_element; i++) {
            HuffmanCode huffman_code = huffman_code_map[((char *) array)[i]];
            bs_write(bit_stream,huffman_code.bit_length,(uint32_t)huffman_code.value);
        }
    } else if (array_element_type == IntArray) {
        for (int i = 0; i < num_element; i++) {
            HuffmanCode huffman_code = huffman_code_map[((int *) array)[i]];
            bs_write(bit_stream,huffman_code.bit_length,(uint32_t)huffman_code.value);
        }
    }	
}
static void codeHuffmanDictionary(struct HuffmanTree *huffman_tree, bs_t *out_stream, const int leaf_node_code_length){
    codeHuffmanNodeBinary(huffman_tree,huffman_tree->root_pos, out_stream, leaf_node_code_length);    
}

static void codeHuffmanNodeBinary(struct HuffmanTree *huffman_tree,size_t node, bs_t *bit_stream, const int leaf_node_code_length) {
    if (huffman_tree->iNode[node].left == -1 && huffman_tree->iNode[node].right == -1) {
        bs_write1(bit_stream,1);
        bs_write(bit_stream,leaf_node_code_length,(uint32_t)huffman_tree->iNode[node].value);

    } else {
        bs_write1(bit_stream,0);
        if (huffman_tree->iNode[node].left != -1) {
            codeHuffmanNodeBinary(huffman_tree, huffman_tree->iNode[node].left, bit_stream, leaf_node_code_length);
        }
        if (huffman_tree->iNode[node].right != -1) {
            codeHuffmanNodeBinary(huffman_tree, huffman_tree->iNode[node].right, bit_stream, leaf_node_code_length);
        }

    }
}

int encodeWithHuffman(const short *data, const int data_length, const int max_symbol_category,
                      HuffmanEncodeContext *huffman_encode_context, bs_t *bitstream){
	int huffman_leaf_length;
    huffman_leaf_length=ceil(log2(max_symbol_category));

    bs_align_1(bitstream);
    bs_write(bitstream, ByteLength * 3,(uint32_t)data_length);
    if (0==data_length){
        return 1;
    }

    int * frequencies = (int *)malloc(max_symbol_category * sizeof(int));
    if (!frequencies)
        return -1;
    for (int i = 0; i < max_symbol_category; i++) {
        frequencies[i] = 0;
    }

    generateShortFrequencyMap(data, data_length, max_symbol_category, frequencies);

    int ret = getHuffmanCodeDict(huffman_encode_context, frequencies, max_symbol_category);

    if (ret==1){
		codeHuffmanDictionary(huffman_encode_context->huffman_tree, bitstream, huffman_leaf_length);		
        codeHuffmanBitStream(data, data_length, ShortArray, huffman_encode_context->huffman_code, bitstream);
        free(frequencies);
        return 1;
    }
    else {
        free(frequencies);
        return -1;
    }
}

static void generateShortFrequencyMap(const short *data, const int data_length, const int max_symbol_category, int *frequency) {
    for (int i = 0; i < data_length; i++) {
        if (data[i]>=0 && data[i]<max_symbol_category){
            frequency[data[i]]++;
        }
    }
}


