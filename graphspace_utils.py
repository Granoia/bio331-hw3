#!/usr/bin/python

## utility functions for posting graphs to graphspace
## BIO331
## Anna Ritz
import sys
import json
import subprocess
__docformat__ = 'reStructuredText'

## URL for original GraphSpace
##URL="http://graphspace.org"
## URL for Reed GraphSpace
URL="http://ec2-52-41-252-78.us-west-2.compute.amazonaws.com/"

def postGraph(graphid,jsonfile,user,password,logfile=None):
	"""
	Posts a graph in 'jsonfile' with id 'graphid' to the account of the user 'user' to GraphSpace.

	:param graphid: string -- ID of GraphSpace graph
	:param jsonfile: string -- JSON file to post. 
	:param user: string -- graph owner's username
	:param password: string -- graph owner's password
	:param logfile: filename for command outputs.  Optional.
	"""
	if logfile:
		logout = open(logfile,'w')
	else:
		logout = None
		
	# check to see if this graph is already in GraphSpace.
	graph_exists = False
	cmd = _constructExistsCommand(graphid,user,password)
	outstring = execute(cmd,logout)
	outstring = json.loads(outstring)
	if outstring["StatusCode"] == 200:
		# a status code of 200 indicates that a graph already exists.
		graph_exists = True

	if graph_exists:  
		#print('\nUpdating existing graph %s from user %s' % (graphid,user))
		cmd = _constructUpdateCommand(graphid,user,password,jsonfile)
		execute(cmd,logout)
	else:
		#print('\nGraph does not exist. Posting new graph %s from user %s' % (graphid,user))
		cmd = _constructPostCommand(graphid,user,password,jsonfile)
		execute(cmd,logout)
	
	if logout:
		print('command output written to %s' % (logfile))
		logout.close()
	return

def deleteGraph(graphid,user,password):
	"""
	Removes a graph (denoted by graphid and user) from GraphSpace.

	:param graphid: ID of GraphSpace graph 
	:param user: graph owner's username
	:param password: graph owner's password
	"""
	print('\nRemoving existing graph %s from user %s' % (graphid,user))
	cmd = _constructDeleteCommand(graphid,user,password)
	execute(cmd)
   
def shareGraph(graphid,user,password,group,group_owner):
	"""
	Shares an existing graph with a group.

	:param graphid: ID of GraphSpace graph
	:param user: graph owner's username
	:param password: graph owner's password
	:param group: group to share graph with.
	:param group_user: owner of group. 
	"""
	print('\nSharing existing graph %s from user %s with group %s owned by %s' % (graphid,user,group,group_owner))
	cmd = _constructShareCommand(graphid,user,password,group,group_owner)
	outstring = execute(cmd)

def unShareGraph(graphid,user,password,group,group_owner):
	"""
	Un-Shares a graph with a group.  Graph is not deleted, but others
	in the group may no longer view it.

	:param graphid: ID of GraphSpace graph
	:param user: graph owner's username
	:param password: graph owner's password
	:param group: group to share graph with.
	:param group_user: owner of group. 
	"""
	print('\nUn-Sharing existing graph %s from user %s with group %s owned by %s' % (graphid,user,group,group_owner))
	cmd = _constructUnShareCommand(graphid,user,password,group,group_owner)
	outstring = execute(cmd)

def makeGraphPublic(graphid,user,password):
	"""
	Makes a graph publicly viewable.

	:param graphid: ID of GraphSpace graph 
	:param user: graph owner's username
	:param password: graph owner's password
	"""
	print('\nMaking graph %s from user %s public' % (graphid,user))
	cmd = _constructPublicGraphCommand(graphid,user,password)
	outstring = execute(cmd)


def makeGraphPrivate(graphid,user,password):
	"""
	Makes a graph privately viewable (makes it no longer public).

	:param graphid: ID of GraphSpace graph 
	:param user: graph owner's username
	:param password: graph owner's password
	"""
	print('\nMaking graph %s from user %s private' % (graphid,user))
	cmd = _constructPrivateGraphCommand(graphid,user,password)
	outstring = execute(cmd)

def makeGraphsWithTagPublic(user,password,tag):
	"""
	Makes all graphs with a tag publicly viewable.

	:param user: graph owner's username
	:param password: graph owner's password
	:param tag: graph tag.
	"""
	print('\nMaking graphs with tag %s from user %s public' % (tag,user))
	cmd = _constructPublicTagCommand(user,password,tag)
	outstring = execute(cmd)  

def makeGraphsWithTagPrivate(user,password,tag):
	"""
	Makes all graphs with a tag privately viewable (make them no longer public).

	:param user: graph owner's username
	:param password: graph owner's password
	:param tag: graph tag.
	"""
	print('\nMaking graphs with tag %s from user %s public' % (tag,user))
	cmd = _constructPrivateTagCommand(user,password,tag)
	outstring = execute(cmd)  

####################################################################
### EXECUTE COMMAND  ###############################################

def execute(cmd,logout=None):
	"""
	Executes the command, using subprocess.Popen.  We need to capture
	the output, so we cannot simply us os.system().

	:param cmd: string -- command to execute
	:param logout: File object -- File of log output or None.
	:return out: string -- output of command
	"""
	print('COMMAND:')
	print(cmd)
	proc = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	(out, err) = proc.communicate()
	## in python3, out is now a byte stream instead of unicde.
	## explicitly cast it. To check if we're using python 2 or 3, see 
	## http://sweetme.at/2013/10/21/how-to-detect-python-2-vs-3-in-your-python-script/
	if sys.version_info > (3,0):
		out =  str( out, encoding='utf8' )
	ind = [i for i in range(len(out)) if out[i]=='}']
	origout = out
	out = out[:ind[0]+1]
	if out != origout:
		print('WARNING: additional HTML detected from output string and trimmed.')
	print('OUTPUT:')
	if 'DOCTYPE' in out or 'html' in out: # HTML is spit out; suppress this and write error statement instead.
		print(out)
		print('ERROR: JSON is improperly formatted.  Check your nodes, edges, or attributes.')
		sys.exit()
	else:
		print(out)
	print('DONE')
	print('')
	if logout:
		logout.write(cmd+'\n')
		logout.write(out+'\n')
	return out

####################################################################
### CURL COMMANDS  #################################################

def _constructExistsCommand(graphid,user,password):
	"""
	Construct curl command to check whether a graph exists.

	:param graphid: string -- ID of GraphSpace graph
	:param user: string -- graph owner username
	:param password: string -- graph owner password
	"""
	cmd = 'curl -X POST %s/api/users/%s/graph/exists/%s/ -F username=%s -F password=%s ; echo' % \
		  (URL,user,graphid,user,password)
	return cmd

def _constructPostCommand(graphid,user,password,jsonfile):
	"""
	Construct curl command to post a graph.

	:param graphid: string -- ID of GraphSpace graph
	:param user: string -- graph owner username
	:param password: string -- graph owner password
	:param jsonfile: string -- JSON file of graph
	"""
	cmd = 'curl -X POST %s/api/users/%s/graph/add/%s/ -F username=%s -F password=%s -F graphname=@%s ; echo'  % \
	   (URL, user, graphid, user, password, jsonfile) 
	return cmd

def _constructUpdateCommand(graphid,user,password,jsonfile):
	"""
	Construct curl command to update (overwrite) a graph.

	:param graphid: string -- ID of GraphSpace graph
	:param user: string -- graph owner username
	:param password: string -- graph owner password
	:param jsonfile: string -- JSON file of graph
	"""
	cmd = 'curl -X POST %s/api/users/%s/graph/update/%s/ -F username=%s -F password=%s -F graphname=@%s ; echo'  % \
		(URL, user, graphid, user, password, jsonfile) 
	return cmd

def _constructDeleteCommand(graphid,user,password):
	"""
	Construct curl command to delete a graph.

	:param graphid: string -- ID of GraphSpace graph
	:param user: string -- graph owner username
	:param password: string -- graph owner password
	"""
	cmd = 'curl -X POST %s/api/users/%s/graph/delete/%s/ -F username=%s -F password=%s ; echo'  % \
		  (URL, user, graphid, user, password) 
	return cmd

def _constructShareCommand(graphid,user,password,groupid,group_owner):
	"""
	Construct curl command to share a graph with a group.

	:param graphid: string -- ID of GraphSpace graph
	:param user: string -- graph owner username
	:param password: string -- graph owner password
	:param groupid: string -- group to share graph with
	:param group_owner: string -- group's owner
	"""
	cmd = 'curl -X POST %s/api/users/graphs/%s/share/%s/%s/ -F username=%s -F password=%s ; echo'  % \
		  (URL, graphid,group_owner,groupid,user,password)
	return cmd

def _constructUnShareCommand(graphid,user,password,groupid,group_owner):
	"""
	Construct curl command to unshare a graph with a group.

	:param graphid: string -- ID of GraphSpace graph
	:param user: string -- graph owner username
	:param password: string -- graph owner password
	:param groupid: string -- group to share graph with
	:param group_owner: string -- group's owner
	"""
	cmd = 'curl -X POST %s/api/users/graphs/%s/unshare/%s/%s/ -F username=%s -F password=%s ; echo'  % \
		  (URL, graphid,group_owner,groupid,user,password)
	return cmd

def _constructPublicGraphCommand(graphid,user,password):
	"""
	Construct curl command to make a graph publicly viewable.

	:param graphid: string -- ID of GraphSpace graph
	:param user: string -- graph owner username
	:param password: string -- graph owner password
	"""
	cmd = 'curl -X POST %s/api/users/%s/graph/makeGraphPublic/%s/ -F username=%s -F password=%s; echo' % \
			(URL,user,graphid,user,password)
	return cmd

def _constructPrivateGraphCommand(graphid,user,password):
	"""
	Construct curl command to make a graph privately viewable.

	:param graphid: string -- ID of GraphSpace graph
	:param user: string -- graph owner username
	:param password: string -- graph owner password
	"""
	cmd = 'curl -X POST %s/api/users/%s/graph/makeGraphPrivate/%s/ -F username=%s -F password=%s; echo' % \
			(URL,user,graphid,user,password)
	return cmd

def _constructPublicTagCommand(user,password,tag):
	"""
	Construct curl command to make all graphs associated with a tag public.

	:param user: string -- graph owner username
	:param password: string -- graph owner password
	:param tag: string -- tag of graphs to make public
	"""
	cmd = 'curl -X POST %s/api/tags/user/%s/%s/makePublic/ -F username=%s -F password=%s; echo' % \
			(URL,user,tag,user,password)
	return cmd

def _constructPrivateTagCommand(user,password,tag):
	"""
	Construct curl command to make all graphs associated with a tag private.

	:param user: string -- graph owner username
	:param password: string -- graph owner password
	:param tag: string -- tag of graphs to make private
	"""
	cmd = 'curl -X POST %s/api/tags/user/%s/%s/makePrivate/ -F username=%s -F password=%s; echo' % \
			(URL,user,tag,user,password)
	return cmd


