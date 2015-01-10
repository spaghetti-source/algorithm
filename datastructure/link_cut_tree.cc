struct link_cut_tree {
  struct node { 
    int x, s; // value and sum
    node *ch[2], *p;
  };
  int sum(node *t) { return t ? t->s : 0; }
  node *update(node *t) {
    if (t) t->s = t->x + sum(t->ch[0]) + sum(t->ch[1]);
    return t;
  }
  node *make_node(int x) { return new node({x, x, 0, 0, 0}); }

  int dir(node *t) { return t != t->p->ch[0]; }
  bool is_root(node *t) {
    return !t->p || (t->p->ch[0] != t && t->p->ch[1] != t);
  }
  void connect(node *p, node *t, int d) {
    p->ch[d] = t; if (t) t->p = p;
  }
  void rot(node *t) {
    node *p = t->p;
    int d = dir(t); 
    if (!is_root(p)) connect(p->p, t, dir(p));
    else             t->p = p->p;
    connect(p, t->ch[!d], d);
    connect(t, p, !d);
    update(p); update(t);
  }
  void splay(node *t) {
    for (; !is_root(t); rot(t))
      if (!is_root(t->p)) rot(dir(t) == dir(t->p) ? t->p : t);
  }
  node *expose(node *t) {
    node *l = 0;
    for (node *s = t; s; s = s->p) {
      splay(s);
      s->ch[1] = l;
      l = update(s);
    }
    splay(t);
    return l;
  }
  void link(node *p, node *t) {
    expose(t);
    t->p = p;
  }
  void cut(node *t) {
    expose(t); 
    t->ch[0] = t->ch[0]->p = 0;
  }
  node *lca(node *s, node *t) {
    expose(s); 
    return expose(t);
  }
  int sum_to_root(node *t) {
    expose(t);
    return sum(t->ch[0]) + t->x;
  }
};
