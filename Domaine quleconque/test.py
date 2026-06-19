l=[]
    s=sum(l)
    while(s!=n):
        l.clear()
        for i in range(0,m):
            l.append(random.randint(0,n))
        print(l)
        print("hello")
        s=sum(l)
    
    print(l)