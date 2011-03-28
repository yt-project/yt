function fill_tree(my_pfs) {
  examine = my_pfs;
  treePanel.root.removeAll();
  Ext.each(my_pfs, function(pf, index) {
      treePanel.root.appendChild(new Ext.tree.TreeNode({text: pf.name,
              leaf:false, expanded:true, iconCls: 'pf_icon'}));
      this_pf = treePanel.root.lastChild
	Ext.each(pf.objects, function(object, obj_index) {
            this_pf.appendChild(new Ext.tree.TreeNode({text: object.name,
		    leaf: true, iconCls: 'data_object'}));
	  });
    });
};
