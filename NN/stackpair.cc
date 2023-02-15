	static struct node
	{ int key; struct node *next; };
	static struct node *heada, *za, *ta;
	static struct node *headb, *zb, *tb;
	stackpairinit() 
	   {
	     heada = (struct node *) malloc(sizeof *heada);
	     za = (struct node *) malloc(sizeof *za);
	     heada->next = za; heada->key=0;
	     za->next = za;
	     za->key = 0;

	     headb = (struct node *) malloc(sizeof *headb);
	     zb = (struct node *) malloc(sizeof *zb);
	     headb->next = zb; headb->key=0;
	     zb->next = zb;
	     zb->key = 0;
	   }
	stackpairflush() 
	   {
             free(heada);
             free(headb);
             free(za);
             free(zb);
	   }
	pushpair(pa, pb)
             int *pa;
             int *pb;
	   {
	     int va;
	     int vb;
	     va = *pa;
	     ta = (struct node *) malloc(sizeof *ta);	
	     ta->key = va; ta->next = heada->next;	
	     heada->next =ta;	

	     vb = *pb;
	     tb = (struct node *) malloc(sizeof *tb);	
	     tb->key = vb; tb->next = headb->next;	
	     headb->next =tb;	
	   }
	poppair(xa,xb)
             int *xa;
             int *xb;
	   {
	     ta = heada->next; heada->next = ta->next;
	     *xa = ta->key;
	     free(ta);
	     tb = headb->next; headb->next = tb->next;
	     *xb = tb->key;
	     free(tb);
	   }
	stackpairempty(i)
            int *i;
	  { 
	    *i = 0;
            if(heada->next == za) *i = 1;
          }
