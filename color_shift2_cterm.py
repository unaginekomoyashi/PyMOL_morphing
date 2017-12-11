from pymol import cmd,stored
import numpy as np
import time,glob

#attention_residues=[324,325,479,502]

def color_shift(start="Apo2",end="cterm2",divide=100):
    stored.start_colors=[]
    cmd.iterate("%s and name ca and chain A" % start,"stored.start_colors.append(color)")
    stored.end_colors=[]
    cmd.iterate("%s and name ca and chain A" % end,"stored.end_colors.append(color)")    
    print stored.start_colors[0],stored.end_colors[0]
    print cmd.get_color_tuple(stored.start_colors[0]),cmd.get_color_tuple(stored.end_colors[0])
    start_color_indice=cmd.get_color_tuple(stored.start_colors[0])
    end_color_indice=cmd.get_color_tuple(stored.end_colors[0])
    print "start ,end" 
    print  start_color_indice,end_color_indice
    print np.array(end_color_indice)-np.array(start_color_indice)
    #print (np.array(end_color_indice)-np.array(start_color_indice))/float(divide)
    graduate_color_indice=(np.array(end_color_indice)-np.array(start_color_indice))/float(divide)
    i=1
    while i<divide:
        i+=1
        #print list(np.array(start_color_indice)+graduate_color_indice*i)
        time.sleep(0.1)
        cmd.refresh()
        cmd.recolor()
        cmd.set_color("temp%s" %i,list(np.array(start_color_indice)+graduate_color_indice*i))
        cmd.color("temp%s" % i,"%s and name ca and chain A" % start)
        
#color_shift(start="cterm",end="cterm2",divide=50)


def paring(start="cterm",end="cterm2",divide=50):
    print "paring : part gradually disaper(start) and appedar(end)"
    cmd.remove("alt B or alt C or alt D")
    cmd.remove("elem H")
    stored.model_chain_resn_resi_name_alt=[]
    cmd.iterate("%s" % (start),"stored.model_chain_resn_resi_name_alt.append((chain,resn,resi,name,alt))")
    #cmd.iterate("%s and byres(n;ca)" % (start),"stored.model_chain_resn_resi_name_alt.append((chain,resn,resi,name,alt))")
    #print stored.model_chain_resn_resi_name_alt
    start_param_list=set(stored.model_chain_resn_resi_name_alt)
    stored.model_chain_resn_resi_name_alt=[]
    cmd.iterate("%s" % (end),"stored.model_chain_resn_resi_name_alt.append((chain,resn,resi,name,alt))")
    #cmd.iterate("%s and byres(n;ca)" % (end),"stored.model_chain_resn_resi_name_alt.append((chain,resn,resi,name,alt))")
    end_param_list=set(stored.model_chain_resn_resi_name_alt)    
    #print start_param_list-end_param_list
    diff_obj_dis=(start_param_list-end_param_list)
    diff_obj_app=(end_param_list-start_param_list)
    #print start_param_list&end_param_list
    common_obj=start_param_list&end_param_list
    
    # start model surplus
    for (chain,resn,resi,name,alt) in diff_obj_dis:
        #print (chain,resn,resi,name,alt)
        cmd.select("diff_dis","vis and chain %s and resi %s and name %s" % (chain,resi,name),merge=1)
    # end model surplus
    for (chain,resn,resi,name,alt) in diff_obj_app:
        #print (chain,resn,resi,name,alt)
        cmd.select("diff_app","vis and chain %s and resi %s and name %s" % (chain,resi,name),merge=1)
    # common part of start and end
    for (chain,resn,resi,name,alt) in common_obj:
        #print (chain,resn,resi,name,alt)
        cmd.select("common_start","%s and chain %s and resi %s and name %s" % (start,chain,resi,name),merge=1)
        cmd.create("common_start","common_start")
    for (chain,resn,resi,name,alt) in common_obj:
        #print (chain,resn,resi,name,alt)
        cmd.select("common_end","%s and chain %s and resi %s and name %s" % (end,chain,resi,name),merge=1)
        cmd.create("common_end","common_end")

paring(start="cterm",end="cterm2",divide=50)

def coord_shift(start="cterm",end="cterm2",divide=10):
    start_model=[(x.chain,x.resi,x.name,x.coord[0],x.coord[1],x.coord[2],x.resn,x.b) for x in cmd.get_model("common_start").atom]
    end_model=[(x.chain,x.resi,x.name,x.coord[0],x.coord[1],x.coord[2],x.resn,x.b) for x in cmd.get_model("common_end").atom]
    ix=0
    while ix<=divide:
        print ix
        fout=open("start%s2_end%s_%s.xyz" % (start,end,ix),"w")
        for i in range(0,len(end_model)):
            print end_model[i]
            chain=end_model[i][0]
            resi=end_model[i][1]
            name=end_model[i][2]
            b=end_model[i][-1]
            resn=end_model[i][-2]
            diffx=float(end_model[i][3]-start_model[i][3])/divide
            diffy=float(end_model[i][4]-start_model[i][4])/divide
            diffz=float(end_model[i][5]-start_model[i][5])/divide
            nowx=float(start_model[i][3])+diffx*ix
            nowy=float(start_model[i][4])+diffy*ix
            nowz=float(start_model[i][5])+diffz*ix
            print nowx,nowy,nowz
            fout.write("%.2f,%.2f,%.2f,%s,%s,%s,%s\n" % (nowx,nowy,nowz,chain,resi,name,resn))
            #if name=="CA" or name=="N" or name=="C" or name=="O":
            #    fout.write("%.2f,%.2f,%.2f,%s,%s,%s,%s\n" % (nowx,nowy,nowz,chain,resi,name,resn))
        fout.close()
        ix+=1
        fout.close()
coord_shift(start="cterm",end="cterm2",divide=50)


def xyz2pdb(start="cterm",end="cterm2",divide=50):
    xyzes=glob.glob("start%s2_end%s_*.xyz" % (start,end))    
    for xyz in xyzes:
        print xyz
        fin=open("%s" % xyz,"r")
        fout=open("%s.pdb" % xyz.split(".")[0],"w")
        ii=1
        for l in fin:
            (x,y,z,c,i,n,r)=l.strip().split(",")
            elem=n[0]
            fout.write("ATOM%7s%5s %s %s%4s    %8.3f%8.3f%8.3f  1.00 %5.2f %s\n" %\
                       (ii,n,r,c,i,float(x),float(y),float(z),float(1),"%s".rjust(12) % elem)) 
            ii+=1
        fin.close()
        fout.close()
xyz2pdb(start="cterm",end="cterm2",divide=50)


def morphing(start="cterm",end="cterm2",divide=50):
    #cmd.show("spheres","i;502,479,325")
    cmd.show("sticks","i;502,479,325,324 and not (n;n,c,o,ca)")
    cmd.show("cartoon")
    pdbs=glob.glob("start%s2_end%s_*.pdb" % (start,end))
    for pdb in pdbs:
        cmd.load(pdb)

morphing()

def color_shift2(start="cterm",end="cterm2",divide=100):
    stored.start_colors=[]
    cmd.iterate("%s and name ca and chain A" % start,"stored.start_colors.append(color)")
    stored.end_colors=[]
    cmd.iterate("%s and name ca and chain A" % end,"stored.end_colors.append(color)")    
    start_color_indice=cmd.get_color_tuple(stored.start_colors[0])
    end_color_indice=cmd.get_color_tuple(stored.end_colors[0])
    graduate_color_indice=(np.array(end_color_indice)-np.array(start_color_indice))/float(divide)
    i=0
    while i<=divide:
        cmd.refresh()
        cmd.recolor()
        cmd.set_color("temp%s" %i,list(np.array(start_color_indice)+graduate_color_indice*i))
        #cmd.color("temp%s" % i,"%s and name ca and chain A and not (e;N,O)" % start)
        cmd.color("temp%s" % i,"start%s2_end%s_%s and chain A and not (e;N,O)" % (start,end,i))
        #cmd.set("sphere_color","temp%s" % i,"start%s2_end%s_%s and  elem C" % (start,end,i))
        cmd.set("stick_color","red","start%s2_end%s_%s and elem O" % (start,end,i))  
        cmd.set("stick_color","blue","start%s2_end%s_%s and elem N" % (start,end,i))         
        cmd.set("cartoon_color","red","start%s2_end%s_%s and chain C" % (start,end,i))
        cmd.set("cartoon_color","tv_blue","start%s2_end%s_%s and chain B" % (start,end,i))
        i+=1
color_shift2(start="cterm",end="cterm2",divide=50)

def morf2(start="cterm",end="cterm2",divide=50):
    objects=[x for x in cmd.get_names("objects") if x.find("start%s2_end%s" % (start,end))==0]
    for obj in objects:
        number=int(obj.split("_")[-1])
        cmd.create("merge_%s_%s" % (start,end),obj,1,number)
        cmd.disable(obj)
    cmd.set("movie_fps",5)
    cmd.show("cartoon")
    cmd.hide("lines")
    #cmd.show("spheres","(i;502,479,325) and not (n:n,c,o,ca)")
    cmd.show("sticks","(i;502,479,325,324)")
    cmd.hide("sticks","(n;n,o)")
    cmd.show("cartoon")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_highlight_color", "white")
morf2()


def setting():
    #cmd.set("bg_color","white")
    cmd.set("cartoon_color","red","chain C")
    cmd.set("cartoon_color","tv_blue","chain B")
    cmd.hide("sticks","(n;n,o)")
    cmd.show("sticks","(i;502,479,325,324)")
    cmd.hide("sticks","(n;n,o,c)")
    util.cbac("resi 9000")
    cmd.show("sticks","resi 9000")
    cmd.set("ray_shadow","off")
setting()



def trans():
    cmd.disable("all")
    objs=[x for x in cmd.get_names("objects") if x.find("start")==0]
    trans=float(1)/len(objs)
    t=0
    for x in objs:
        t+=trans
        cmd.enable(x)
        print x,t
        cmd.set("cartoon_transparency",t)
        cmd.set("stick_transparency",t)       
        cmd.png("%s.png" %x,ray=1)
        cmd.disable(x)
        
trans()





