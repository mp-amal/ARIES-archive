
import json
from datetime import datetime,date,timedelta
import os
import paramiko
import yaml

import json
import smtplib
from datetime import datetime,date
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from datetime import datetime


with open("/home/archive/Documents/ARIES-archive/config/paths.yaml", "r") as file:
    paths_yaml = yaml.safe_load(file)



template = {
    'total': '-',
    'science': '-',
    'flat': '-',
    'bias': '-',
    'lamp': '-',
    'status': 'not completed'
}

def generate_email_html(all_data_for_date, your_name="Amal MP"):
    """
    all_data_for_date: dict loaded from JSON (same structure you showed).
    Returns: HTML string ready to send as an email body.
    """

    data = all_data_for_date

    # -------- 1) Flatten into rows --------
    # Each row: (instrument, sub_instrument, total, science, flat, bias, lamp, status, comment)
    rows = []

    for inst_name, inst_data in data.items():
        if not isinstance(inst_data, dict):
            continue

        inst_status = inst_data.get("status", "-")
        inst_comment = inst_data.get("comment", "")

        # Case A: top-level counts (like ST)
        has_top_counts = any(k in inst_data for k in ("total", "science", "flat", "bias", "lamp"))
        if has_top_counts:
            rows.append((
                inst_name,
                "-",  # no sub-instrument
                inst_data.get("total", "-"),
                inst_data.get("science", "-"),
                inst_data.get("flat", "-"),
                inst_data.get("bias", "-"),
                inst_data.get("lamp", "-"),
                inst_status,
                inst_comment
            ))

        # Case B: nested sub-instruments (like DFOT["2K_IMG1"], DOT["ADFOSC"])
        for sub_name, sub_data in inst_data.items():
            if not isinstance(sub_data, dict):
                continue
            if not any(k in sub_data for k in ("total", "science", "flat", "bias", "lamp")):
                continue

            sub_status = sub_data.get("status", inst_status)
            sub_comment = sub_data.get("comment", inst_comment)

            rows.append((
                inst_name,
                sub_name,
                sub_data.get("total", "-"),
                sub_data.get("science", "-"),
                sub_data.get("flat", "-"),
                sub_data.get("bias", "-"),
                sub_data.get("lamp", "-"),
                sub_status,
                sub_comment
            ))

    # -------- 2) Summary statistics --------
    total_instruments = len(data)

    completed_count = sum(
        1 for v in data.values()
        if isinstance(v, dict) and str(v.get("status", "")).lower() == "completed"
    )
    pending_count = total_instruments - completed_count

    def to_int_or_zero(x):
        try:
            return int(x)
        except Exception:
            return 0

    total_science_frames = sum(to_int_or_zero(r[3]) for r in rows)  # index 3 is "science"

    # -------- 3) Build table rows HTML --------
    def badge_style(status):
        status_str = str(status).lower()
        if status_str == "completed":
            return (
                "display:inline-block; padding:2px 6px; border-radius:999px; "
                "background:#e0fbea; color:#166534; border:1px solid #a7f3d0; "
                "font-size:11px; font-weight:bold; text-transform:uppercase;"
            )
        else:
            return (
                "display:inline-block; padding:2px 6px; border-radius:999px; "
                "background:#fef2f2; color:#b91c1c; border:1px solid #fecaca; "
                "font-size:11px; font-weight:bold; text-transform:uppercase;"
            )

    def esc(x):
        return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    row_html_parts = []
    for (inst, subinst, total, sci, flat, bias, lamp, status, comment) in rows:
        row_html_parts.append(f"""
        <tr>
          <td style="padding:7px 8px; border-bottom:1px solid #e5e7eb;">{esc(inst)}</td>
          <td style="padding:7px 8px; border-bottom:1px solid #e5e7eb;">{esc(subinst)}</td>
          <td style="padding:7px 8px; border-bottom:1px solid #e5e7eb;">{esc(total)}</td>
          <td style="padding:7px 8px; border-bottom:1px solid #e5e7eb;">{esc(sci)}</td>
          <td style="padding:7px 8px; border-bottom:1px solid #e5e7eb;">{esc(flat)}</td>
          <td style="padding:7px 8px; border-bottom:1px solid #e5e7eb;">{esc(bias)}</td>
          <td style="padding:7px 8px; border-bottom:1px solid #e5e7eb;">{esc(lamp)}</td>
          <td style="padding:7px 8px; border-bottom:1px solid #e5e7eb;">
            <span style="{badge_style(status)}">{esc(status)}</span>
          </td>
          <td style="padding:7px 8px; border-bottom:1px solid #e5e7eb;">{esc(comment)}</td>
        </tr>
        """)
    if not row_html_parts:
        row_html_parts.append(
            '<tr><td colspan="9" style="padding:8px; text-align:center;">No data</td></tr>'
        )
    rows_html = "\n".join(row_html_parts)
    inst = 'TANSPEC'

    # -------- 4) Plug into email-safe HTML template --------
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <title>Observation Data Status</title>
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
</head>
<body style="margin:0; padding:0; background-color:#f3f4f6; font-family:Arial, Helvetica, sans-serif; color:#111827;">
  <center style="width:100%; background-color:#f3f4f6; padding:16px 0;">
    <table role="presentation" cellpadding="0" cellspacing="0" width="100%" style="max-width:760px; background:#ffffff; border-radius:12px; border:1px solid #d1d5db;">
      <!-- Header -->
      <tr>
        <td style="padding:18px 20px; background:#1d4ed8; background:linear-gradient(135deg,#1d4ed8,#0f172a); color:#e5e7eb; border-radius:12px 12px 0 0;">
          <h1 style="margin:0; font-size:20px; font-weight:bold;">Observation Data Status Summary ({str(datetime.strptime(yesterday_date, "%Y%m%d").strftime("%Y-%m-%d"))})</h1>
          <p style="margin:6px 0 0 0; font-size:13px; color:#c7d2fe;">
            DFOT • ST • DOT  &nbsp;|&nbsp; Auto-generated daily report
          </p>
        </td>
      </tr>

      <!-- Intro -->
      <tr>
        <td style="padding:16px 20px 4px 20px; font-size:14px; line-height:1.5;">
          <p style="margin:0 0 10px 0;">Dear Team,</p>
          <p style="margin:0 0 10px 0;">
            Please find below the latest consolidated status of the observational data processed
            for DFOT, ST and DOT.
          </p>
        </td>
      </tr>

      <!-- Simple stats row -->
      <tr>
        <td style="padding:6px 20px 14px 20px;">
          <table role="presentation" cellpadding="0" cellspacing="0" width="100%" style="border-collapse:collapse;">
            <tr>
              <td style="padding:8px; font-size:12px; border:1px solid #e5e7eb; border-radius:8px; background:#f9fafb;">
                <strong>Telescope Archiving:</strong> {total_instruments} &nbsp;&nbsp;
                <strong>Completed:</strong> {completed_count} &nbsp;&nbsp;
                <strong>Pending:</strong> {pending_count} &nbsp;&nbsp;
                <strong>Total Science Frames:</strong> {total_science_frames}
              </td>
            </tr>
          </table>
        </td>
      </tr>

      <!-- Main table -->
      <tr>
        <td style="padding:4px 20px 16px 20px;">
          <table cellpadding="0" cellspacing="0" width="100%" style="border-collapse:collapse; border:1px solid #e5e7eb; border-radius:10px; overflow:hidden; font-size:12px;">
            <tr style="background-color:#eff2ff;">
              <th align="left" style="padding:8px 8px; border-bottom:1px solid #e5e7eb;">Telescope</th>
              <th align="left" style="padding:8px 8px; border-bottom:1px solid #e5e7eb;">Instrument</th>
              <th align="left" style="padding:8px 8px; border-bottom:1px solid #e5e7eb;">Total</th>
              <th align="left" style="padding:8px 8px; border-bottom:1px solid #e5e7eb;">Science</th>
              <th align="left" style="padding:8px 8px; border-bottom:1px solid #e5e7eb;">Flat</th>
              <th align="left" style="padding:8px 8px; border-bottom:1px solid #e5e7eb;">Bias</th>
              <th align="left" style="padding:8px 8px; border-bottom:1px solid #e5e7eb;">Lamp</th>
              <th align="left" style="padding:8px 8px; border-bottom:1px solid #e5e7eb;">Status</th>
              <th align="left" style="padding:8px 8px; border-bottom:1px solid #e5e7eb;">Comment</th>
            </tr>
            {rows_html}
          </table>
        </td>
      </tr>

      <!-- Closing -->
      <tr>
        <td style="padding:10px 20px 6px 20px; font-size:13px; line-height:1.5;">
          <p style="margin:0 0 8px 0;">
            Not completed, corresponding data folder is not available in the telescope-data user server for processing.
          </p>

        </td>
      </tr>

      <!-- Footer -->
      <tr>
        <td style="padding:8px 20px 12px 20px; font-size:11px; color:#6b7280; background:#f3f4f6; border-radius:0 0 12px 12px;">
          This summary is auto-generated from the latest JSON object for the last night.
        </td>
      </tr>
    </table>
  </center>
</body>
</html>
"""
    return html

files = [paths_yaml['DFOT_JSON'], paths_yaml['ST_JSON'],paths_yaml['DOT_JSON']]

today_date = date.today().strftime("%Y%m%d")
yesterday_date = (date.today() - timedelta(days=1)).strftime("%Y%m%d")
date_key = str(yesterday_date)
year = today_date[:4]
# date_key ='20250228'
all_data_for_date = {}

for fname in files:
    with open(fname, "r") as f:
        data = json.load(f)

    # store data for this file under its name (without .json)
    base = fname.rsplit("_", 1)[0]
    telescope = os.path.basename(base)
    # print(telescope)
    if date_key in data:
        all_data_for_date[telescope] = data[date_key]
    else:
        all_data_for_date[telescope] = None  # or skip / handle as you like

# print(f"Data for date {date_key} from all files:")


for teles in ['ST','DFOT','DOT']:
    if all_data_for_date.get(teles) is None:

        all_data_for_date[teles] = template.copy()
        all_data_for_date[teles]['status'] = 'No observation'
        all_data_for_date[teles]['comment'] = 'No data'

    else:
        all_data_for_date[teles]['status'] = 'completed'
        all_data_for_date[teles]['comment'] = '✅'

# print(all_data_for_date)
html_body = generate_email_html(all_data_for_date, your_name="Amal MP")
# print(html_body)

def send_email(today_date,subject, to_email, html_body,cc_emails):
    # Email settings
    sender_email = "ariesdataarchive@gmail.com"
    sender_password = "ccxp ixke fgtb bhmt"
    smtp_server = "smtp.gmail.com"
    smtp_port = 587

  
    # Create a multipart email
    msg = MIMEMultipart("alternative")
    msg['From'] = sender_email
    msg['To'] = to_email
    msg['Subject'] = subject
    if cc_emails:
        # If multiple CCs are provided, join them into a single string
        msg['Cc'] = ', '.join(cc_emails)  # Add multiple CCs
    # Attach plain text body (optional)
    # if plain_body:
    #     msg.attach(MIMEText(plain_body, 'plain'))

    # Attach the dynamic HTML body
    msg.attach(MIMEText(html_body, 'html'))

    # Combine To and Cc for sending
    recipients = [to_email]
    if cc_emails:
        recipients += cc_emails  # Add all CCs to the recipients list
    print(recipients)
    try:
        # Connect to the SMTP server and send the email
        server = smtplib.SMTP(smtp_server, smtp_port)
        server.starttls()  # Secure the connection
        server.login(sender_email, sender_password)
        text = msg.as_string()
        server.sendmail(sender_email, recipients, text)
        server.quit()

        print("Email sent successfully!")
    except Exception as e:
        print(f"Failed to send email: {str(e)}")

today_date = date.today().strftime("%Y%m%d")
year = today_date[:4]
subject = "Data Archive Daily report : "+date.today().strftime("%Y/%m/%d")  
# cc_emails=['thisisamalmp@gmail.com','sanjit@aries.res.in','neelam@aries.res.in','govind@aries.res.in','suvendu@aries.res.in','ashishdevaraj@aries.res.in','geetarangwal@aries.res.in']
cc_emails =['thisisamalmp@gmail.com']
# to_email = "kuntal@aries.res.in"
to_email= 'amalmpdiaries@gmail.com'
# to_email ='thisisamalmp@gmail.com'
# cc_emails =['amalmpdiaries@gmail.com']
send_email(today_date,subject, to_email, html_body,cc_emails)
