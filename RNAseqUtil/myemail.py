#!/usr/bin/python

def email2me(subject,body,attachment=[]):
	src = 'luzhao1986@gmail.com'
	rec = [src]
	sendEmail(src,rec,subject,body,attachment)

def sendEmail(src,receivers,subject,body,attachment = [],user = 'luzhao1986',passwd = 'cheniong'):
	
	import os
	import smtplib
	from email.mime.text import MIMEText
	from email.MIMEMultipart import MIMEMultipart
	from email.Utils import formatdate
	from email.MIMEBase import MIMEBase
	from email import Encoders
	import gzip
	msg = MIMEMultipart()
	msg.attach(MIMEText(body)) 
	to = receivers 
	msg['Subject'] = subject
	msg['From'] = src 
	msg['To'] = ','.join(receivers) 
	msg['Date'] = formatdate(localtime=True)
	if attachment:
		list = ' '.join(attachment)
		os.system('tar czvf attachment.tar.gz ' + list)
		part = MIMEBase('application','zip')
		part.set_payload(gzip.open('attachment.tar.gz','r').read())
		Encoders.encode_base64(part)
		part.add_header('Content-Disposition','attachment', filename = os.path.basename('attachment.tar.gz'))
		msg.attach(part)	
	try:
		s = smtplib.SMTP('smtp.gmail.com:587')
		s.ehlo()
		s.starttls()
		s.login(user,passwd)
		s.sendmail(src,to,msg.as_string())
		s.close()
	except Exception as error:
		print error
		print "Error: unable to send email"
def parse_email(file):
	content = {}
	for line in file:
		rec = line.rstrip('\n').split(':')
		if rec[0] == 'subject':
			content['subject'] = rec[1]
		elif rec[0] == 'body':
			content['body']=rec[1]
		elif rec[0] == 'attachment':
			content['attachment'] = rec[1:]
	return content	
if __name__ == '__main__':
	import sys
	import os
	content = parse_email(sys.stdin)
	src = 'luzhao1986@gmail.com'
	receivers = [src]
	sendEmail(src,receivers,content['subject'],content['body'],content['attachment'])
