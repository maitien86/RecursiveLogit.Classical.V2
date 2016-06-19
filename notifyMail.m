function notifyMail( mode , varargin)
% Email notification script for MATLAB
% Usage: To sent an email from any address to any address via Matlab
%        script. Can be used to notify the state of the running script to
%        the programmer's inbox. 
%
% To initialise, run the command:
%    notifyMail('set',address,password,server,port);
%        address  - sender's email account address 
%        password - sender's password
%        server   - smtp server, use gmail server by default
%        port     - smtp port, use gmail port by default
%
% To send email, run the command:
%    notifyMe('send',text);
%        text     - a body of the email

global myaddress;

if strcmp(mode,'set')    
    port = '465';
    smtp_server = 'smtp.gmail.com';
    
    if (size(varargin) == size([1 2 3 4]))
        smtp_server = varargin(1,3);
        port = varargin(1,4);
    end
    
    if (size(varargin) == size([1 2 3]))
        smtp_server = varargin(1,3);
    end
        
    if(size(varargin,2) >= 2)
        myaddress = varargin(1,1);
        mypassword = varargin(1,2);
        
        setpref('Internet','E_mail',myaddress);
        setpref('Internet','SMTP_Server',smtp_server);
        setpref('Internet','SMTP_Username',myaddress);
        setpref('Internet','SMTP_Password',mypassword);

        props = java.lang.System.getProperties;
        props.setProperty('mail.smtp.auth','true');
        props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
        props.setProperty('mail.smtp.socketFactory.port',port);
    else
        disp('Email Username and password must be defined');
    end
    
    
elseif strcmp(mode,'send')
    
    subject = 'Notification mail from MATLAB';
    text = ['This email is a notification sent from ', myaddress 10 ...
        'Sent from Matlab'];

%     if (size(varargin) == [1,2])
%         subject = varargin(1,1);
%         text = [varargin(1,2) 10 'Sent from Matlab'];
%     end

    if (size(varargin) == [1,1])
        
        text = [varargin(1,1) 10 'Sent from Matlab'];
    end

    sendmail(myaddress, subject, text);
    
else
    disp('Unrecognised mode');
end
