/* Copyright (c) 2011 the authors listed at the following URL, and/or
the authors of referenced articles or incorporated external code:
http://en.literateprograms.org/Disjoint_set_data_structure_(C)?action=history&offset=20080516180553

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Retrieved from: http://en.literateprograms.org/Disjoint_set_data_structure_(C)?oldid=13366
*/

#include <stdlib.h>

#include "union_find.h"

forest_node* MakeSet(void* value) {
    forest_node* node = malloc(sizeof(forest_node));
    node->value = value;
    node->parent = NULL;
    node->rank = 0;
    return node;
}

void Union(forest_node* node1, forest_node* node2) {
    if (node1->rank > node2->rank) {
        node2->parent = node1;
    } else if (node2->rank > node1->rank) {
        node1->parent = node2;
    } else { /* they are equal */
        node2->parent = node1;
        node1->rank++;
    }
}

forest_node* Find(forest_node* node) {
    forest_node* temp;
    /* Find the root */
    forest_node* root = node;
    while (root->parent != NULL) {
        root = root->parent;
    }
    /* Update the parent pointers */
    while (node->parent != NULL) {
        temp = node->parent;
        node->parent = root;
        node = temp;
    }
    return root;
}


