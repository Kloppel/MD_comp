def loading_bar(current: float, total: float, msg="",
                bar_length=20, bar_token="=",
                empty_token=" ", pointer_token=">",
                finish_token="X"):
    progress = current/total
    bar_length-=1
    if progress >= 1:
        pointer_token = finish_token
    bar = f"[{bar_token*int(progress*bar_length)}{pointer_token}{empty_token*(bar_length-int(progress*bar_length))}] {progress*100:.2f}% Complete"
    return f"\r{msg}{bar}"

def seconds_to_time(seconds: float):
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return hours, minutes, seconds

def seconds_to_timestring(seconds: float):
    hours, minutes, seconds = seconds_to_time(seconds)
    time_string=""
    if hours > 0:
        time_string+=f"{hours} hours, "
    if minutes > 0:
        time_string+=f"{minutes} minutes, "
    time_string+=f"{seconds} seconds."
    return time_string

def get_keyString(d: dict):
    keystring=""
    for key in d.keys():
        keystring+=f"{key}, "
    return keystring[:-2]


