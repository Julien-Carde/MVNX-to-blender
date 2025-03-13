# MVNX-to-blender
MVNX Motion Capture to Blender Animation
This script processes MVNX motion capture files, extracts joint angles, and animates a Blender armature by applying motion data frame by frame. It served as an initial experiment for a project aimed at developing a web application that allows patients to view a real-time 3D avatar of themselves, generated from an experimental set of sensors placed on their body.

Features
	•	Parses MVNX motion capture files (.mvnx).
	•	Extracts joint angles, positions, and velocities from recorded motion.
	•	Detects clapping events using signal processing.
	•	Maps MVNX joints to Blender armature bones.
	•	Animates a Blender rig using motion capture data.

 The script extracts joint angles and applies them to Blender bones:
	•	It reads jointAngle data from the MVNX file.
	•	Converts the angles to Blender-compatible Euler rotations.
	•	Inserts keyframes for each frame to create the animation.
