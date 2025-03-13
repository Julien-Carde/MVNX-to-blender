import untangle
import matplotlib
import matplotlib.pyplot as plt
import scipy.signal
import copy
import numpy as np

def plotData(data,tag,labels=None, dataRange=None, secondsRange=None,
             showDataPlot=True):
    total=len(data)
    info=[]
    
    for i in range(total): info.append(data[i][tag])
    
    if labels==None:
        labels=list(map(lambda x: str(x), list(range(len(info[0])))))
    
    if secondsRange!=None:
        dataRange=list(map(lambda x:int(x*60), secondsRange))

    if dataRange!=None:
        info=info[dataRange[0]:dataRange[1]]
    else:
        dataRange=[0, len(data)]
        
    separated={}
    for i in range(len(info)):
        for j in range(len(info[0])):
            if i==0:
                separated[labels[j]]=[info[i][j]]
            else:
                separated[labels[j]].append(info[i][j])
    
    t=list(range(dataRange[0],dataRange[1]))
    
    if showDataPlot:
        fig=plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        for i in separated:
            ax.plot(t,separated[i], label=i)
        plt.show()
    
    return separated

def findClaps(data,tag,n_claps=2,labels=None, dataRange=None, secondsRange=None, spikeSep=10,
              showDataPlot=True, drop=None):
    separated=plotData(data,tag,labels=labels, dataRange=dataRange, secondsRange=secondsRange,
                       showDataPlot=showDataPlot)
    #lets found the spikes in each data
    fcs=[]
    for i in separated:
        globalMax=max(separated[i])
        fc=[]
        count=0
        factor=0.01
        while len(fc)!=n_claps and count<10000:
            x=copy.deepcopy(separated[i])
#            
#            x=list(map(lambda x: x if x>globalMax*factor else globalMax*factor,x))
#            x=np.array(x)
#            fc=scipy.signal.argrelextrema(x,np.greater,order=10)
            fc,_=scipy.signal.find_peaks(x,height=globalMax*factor, distance=spikeSep)
            
            factor*=1.05
            if factor>=1:
                break
            count+=1
        
        fcs.append([i,fc])
    
    n=n_claps
    if drop!=None:
        if type(drop)==int:
            fcs.remove(fcs[drop])
        else:
            cpfcs=[]
            for i in range(len(fcs)):
                if i not in drop:
                    cpfcs.append(fcs[i])
            fcs=cpfcs
                
                
    if not all(len(i[1])==n for i in fcs):
        goods=[]
        for j in fcs:
            if len(j[1])==n:
                goods.append(j[1])
        if len(goods)==0:  
            raise ValueError ('no solution')
        else:
            claps=sum([i for i in goods])
            res=[int(i/len(goods)) for i in claps] 
    else:
        claps=sum([i[1] for i in fcs])
        res=[int(i/len(fcs)) for i in claps]
    
    res=[res[i]+dataRange[0] for i in range(len(res))]
    return res


def get_default_segment_index():
    index=['pelvis','l5','l3','t12','t8','neck','head','right_shoulder',
           'right_upper_arm','right_forearm','right_hand','left_shoulder',
           'left_upper_arm','left_forearm','left_hand','right_upper_leg',
           'right_lower_leg','right_foot','right_toe','left_upper_leg',
           'left_lower_leg','left_foot','left_toe']
    return index
def get_default_sensor_index():
    index=['pelvis','t8','head','right_shoulder','right_upper_arm','right_fore_arm',
           'right_hand','left_shoulder','left_upper_arm','left_fore_arm', 'left_hand',
           'right_upper_leg','right_lower_leg','right_foot',
           'left_upper_leg','left_lower_leg','left_foot']
    return index
def get_default_joint_index():
    index=['L5S1','L4L3','L1T12','T9T8','T1C7','C1Head','RightC7Shoulder',
           'RightShoulder','RightElbow','RightWrist','LeftC7Shoulder',
           'LeftShoulder','LeftElbow','LeftWrist', 'RightHip','RightKnee',
           'RightAnkle','RightBallFoot','LeftHip','LeftKnee','LeftAnkle',
           'LeftBallFoot']
    index=list(map(lambda x:x+'j',index))
    return index
def get_dictionary_from_class(obj):
    ll_atrib=dir(obj)
    d={}
    for attribute in ll_atrib:
        d[attribute]=getattr(obj,attribute)
    return d

def prepareline(line, tag, index, scaled):
    linep=line.split(' ')
    # Only the information of sensor, segments or joints pass through this function
    #First two lines: segments, 3rd line: sensors, 
    how_many_info={'position':3, 'orientation':4, 'velocity':3, 'acceleration': 3,
                   'angularVelocity': 3, 'angularAcceleration':3, 
                   'sensorMagneticField':3,'sensorOrientation': 4, 'sensorFreeAcceleration':3, 
                   'jointAngle':3, 'jointAngleXZY':3,'jointAngleErgo':3,
                   'centerOfMass':3, 'footContacts':1}
                    
    linep=list(map(lambda x: scaled*float(x), linep)) #pass to float
    if tag in how_many_info:
        n=how_many_info[tag]
    else:
        s='No tag <'+tag+'> in how_many_info'
        raise ValueError (s)
        
    linep=[linep[x:x+n] for x in range(0, len(linep),n)]
    if len(linep)!=len(index):
        print(index)
        s='Error using '+tag
        raise ValueError (s)
    coord={}
    for i in range(len(index)):coord[index[i]]=linep[i]
    return coord
class MVNX_Index():
    def __init__(self, info, prop=False):
        self.segment_index=get_default_segment_index()
        self.sensor_index=get_default_sensor_index()
        self.joint_index=get_default_joint_index()
        self.footContacts_index=['LeftFoot_Heel', 'LeftFoot_Toe','RightFoot_Heel', 'RightFoot_Toe']
        self.mass_index=['Center_Of_Mass']
        self.ergo_index=['jnt1','jnt2','jnt3','jnt4']
        self.info=info
        if prop:
            for index in [self.segment_index, self.sensor_index, self.joint_index]:
                index.append('prop')
    def __getitem__(self, tag):
        if tag in self.info:
            if tag in ['position','orientation','velocity','acceleration',
                   'angularVelocity', 'angularAcceleration']:
                return self.segment_index
            elif tag in ['sensorMagneticField','sensorOrientation','sensorFreeAcceleration']:
                return self.sensor_index
            elif tag in ['jointAngle', 'jointAngleXZY']:
                return self.joint_index
            elif tag=='jointAngleErgo':
                return self.ergo_index
            elif tag=='footContacts':
                return self.footContacts_index
            elif tag=='centerOfMass':
                return self.mass_index
            else:
                return None
class MVNX():
    def __init__(self,mvnx_root):
        self.root=mvnx_root
        self.xml_file=untangle.parse(mvnx_root)
        self.frames=self.xml_file.mvnx.subject.frames.frame
        self.frames=self.frames[3:]
        self.total_frames=len(self.frames)
        self.available_info=[i for i in get_dictionary_from_class(self.frames[3])]
        self.all_index=MVNX_Index(self.available_info)
        
        self.mvn_version=self.xml_file.mvnx.mvn['version']
        self.mvn_build=self.xml_file.mvnx.mvn['build']
        self.subject=self.xml_file.mvnx.subject._attributes

    def get_info(self,tag, scaled=1):
        if tag not in self.available_info:
            raise ValueError ('Not available info in MVNX')
        else:
            ll=[]
            for i in range(self.total_frames):
                s=getattr(self.frames[i],tag)
                s=getattr(s,'cdata')
                d=prepareline(s,tag, self.all_index[tag], scaled)
                ll.append(d)
        return ll
class MVNX_1_PROP():
    def __init__(self,mvnx_root):
        self.root=mvnx_root
        self.xml_file=untangle.parse(mvnx_root)
        self.frames=self.xml_file.mvnx.subject.frames.frame
        self.frames=self.frames[3:]
        self.total_frames=len(self.frames)
        self.available_info=[i for i in get_dictionary_from_class(self.frames[3])]
        self.all_index=MVNX_Index(self.available_info, prop=True)

        self.mvn_version=self.xml_file.mvnx.mvn['version']
        self.mvn_build=self.xml_file.mvnx.mvn['build']
        self.subject=self.xml_file.mvnx.subject._attributes

    def get_info(self,tag, scaled=1):
        if tag not in self.available_info:
            raise ValueError ('Not available info in MVNX')
        else:
            ll=[]
            for i in range(self.total_frames):
                s=getattr(self.frames[i],tag)
                s=getattr(s,'cdata')
                d=prepareline(s,tag, self.all_index[tag], scaled)
                ll.append(d)
        return ll

def obtainMVNX(mvnx_root):
    mvnx_file=MVNX_1_PROP(mvnx_root)
    if mvnx_file.xml_file.mvnx.subject['segmentCount']=='23':
        mvnx_file=MVNX(mvnx_root)
    return mvnx_file



##############################################################################

import bpy

# New function to animate the armature based on MVNX file data
def animate_armature_from_mvnx(armature_name, mvnx_file_path, bone_joint_mapping):
    # Load MVNX file
    mvnx_data = MVNX(mvnx_file_path)

    # Get joint angle data
    joint_angles = mvnx_data.get_info('jointAngle')

    # Ensure the armature exists
    if armature_name not in bpy.data.objects:
        print(f"Armature '{armature_name}' not found.")
        return

    # Get the armature
    armature = bpy.data.objects[armature_name]

    # Clear existing keyframes for all bones in the armature
    if armature.animation_data and armature.animation_data.action:
        for fcurve in armature.animation_data.action.fcurves:
            armature.animation_data.action.fcurves.remove(fcurve)

    # Iterate over each bone and apply joint data
    for bone_name, joint_name in bone_joint_mapping.items():
        # Skip bones with placeholder or undefined names
        if 'YourBlenderBoneNameFor' in bone_name:
            print(f"Skipping undefined bone mapping for joint '{joint_name}'")
            continue

        if bone_name not in armature.pose.bones:
            print(f"Bone '{bone_name}' not found in armature '{armature_name}'.")
            continue
        bone = armature.pose.bones[bone_name]

        # Iterate over each frame and apply the joint angle
        for frame_number, frame_data in enumerate(joint_angles):
            if joint_name not in frame_data:
                print(f"Joint '{joint_name}' not found in frame {frame_number}: available joints are {frame_data.keys()}")
                continue


            # Retrieve the joint angle for the specific joint
            joint_angle = frame_data[joint_name]

            # Assuming joint_angle is a list or tuple of [X, Y, Z] Euler angles
            x_angle, y_angle, z_angle = joint_angle

            # Set the frame in Blender
            bpy.context.scene.frame_set(frame_number)

            # Set rotation and insert keyframe
            bone.rotation_mode = 'XYZ'
            bone.rotation_euler = (x_angle, y_angle, z_angle)
            bone.keyframe_insert(data_path="rotation_euler")



armature_name = "Armature"
mvnx_file_path = "/Users/juliencardeillac/Downloads/example.mvnx"

bone_joint_mapping = {
    'YourBlenderBoneNameForL5S1': 'L5S1j',
    'YourBlenderBoneNameForL4L3': 'L4L3j',
    'YourBlenderBoneNameForL1T12': 'L1T12j',
    'YourBlenderBoneNameForT9T8': 'T9T8j',
    'YourBlenderBoneNameForT1C7': 'T1C7j',
    'mixamorig:Head': 'C1Headj',
    'YourBlenderBoneNameForRightC7Shoulder': 'RightC7Shoulderj',
    'mixamorig:RightShoulder': 'RightShoulderj',
    'mixamorig:RightForeArm': 'RightElbowj',
    'mixamorig:': 'RightWristj',
    'mixamorig:': 'LeftC7Shoulderj',
    'YourBlenderBoneNameForLeftShoulder': 'LeftShoulderj',
    'YourBlenderBoneNameForLeftElbow': 'LeftElbowj',
    'mixamorig:LeftHand': 'LeftWristj',
    'YourBlenderBoneNameForRightHip': 'RightHipj',
    'YourBlenderBoneNameForRightKnee': 'RightKneej',
    'YourBlenderBoneNameForRightAnkle': 'RightAnklej',
    'YourBlenderBoneNameForRightBallFoot': 'RightBallFootj',
    'YourBlenderBoneNameForLeftHip': 'LeftHipj',
    'YourBlenderBoneNameForLeftKnee': 'LeftKneej',
    'YourBlenderBoneNameForLeftAnkle': 'LeftAnklej',
    'YourBlenderBoneNameForLeftBallFoot': 'LeftBallFootj',
}

#animate_armature_from_mvnx(armature_name, mvnx_file)
animate_armature_from_mvnx(armature_name, mvnx_file_path, bone_joint_mapping)





