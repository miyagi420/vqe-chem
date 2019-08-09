if __name__ == '__main__':
    q_name = 'QVM'
    qvm = pyquil.api.QVMConnection()
    if q_name == 'QVM':
        q_device = qvm
    else:
        q_device = pyquil.api.QPUConnection(q_name)
