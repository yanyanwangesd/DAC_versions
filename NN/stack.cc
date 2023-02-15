	static struct node
	{ int key; struct node *next; };
	static struct node *head, *z, *t;
	stackinit() 
	   {
	     head = (struct node *) malloc(sizeof *head);
	     z = (struct node *) malloc(sizeof *z);
	     head->next = z; head->key=0;
	     z->next = z;
	     z->key = 0;
	   }
	push(p)
           int *p;
	   {
	     int v;
	     v = *p;
	     t = (struct node *) malloc(sizeof *t);	
	     t->key = v; t->next = head->next;	
	     head->next =t;	
	   }
	pop(x)
           int *x;
	   {
	     t = head->next; head->next = t->next;
	     *x = t->key;
	     free(t);
	   }
	stackempty(i)
          int *i;
	  { 
	    *i = 0;
            if(head->next == z) *i = 1;
          }
        stackflush()
           {
             free(head);
             free(z);
           }
